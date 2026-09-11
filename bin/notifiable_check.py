#!/usr/bin/env python3
"""
notifiable_check.py - screen a taxonomic profile against the notifiable-disease organism
list (AU / UK / US) and attach an explicit confidence assessment to every hit.

Why this exists
---------------
On 2026-09-11 an interim report named "Treponema" in a patient without immediately
distinguishing *T. pallidum* (syphilis - notifiable) from the oral treponemes actually
present (*T. denticola*, *T. medium*, *T. vincentii* - not notifiable, and expected in a
periodontal infection). The genus name alone carried an implication the data did not
support.

The trap generalises badly, because the notifiable list contains genera that are dominated
by commensals in respiratory samples: Streptococcus, Neisseria, Haemophilus, Campylobacter,
Clostridium, Escherichia, Klebsiella, Staphylococcus. Flagging on genus alone would fire on
almost every specimen; ignoring genus would have missed the Treponema question entirely.

So this tool separates three cases and never collapses them:

  SPECIES_NOTIFIABLE  the call *is* a notifiable species            -> validate before reporting
  GENUS_UNRESOLVED    genus contains a notifiable species, but the
                      call is unresolved (e.g. "Treponema sp.")     -> MUST resolve first
  GENUS_CONGENER      genus contains a notifiable species, but this
                      call is a different, non-notifiable species   -> record, do not alarm

and it refuses to clear anything it has not seen evidence for.

CALIBRATION WARNING
-------------------
The confidence thresholds below are NOT calibrated. They were set by eye from a single
cohort (7 ONT sputum metagenomes) and are not derived from a labelled truth set. There is
no sensitivity/specificity estimate behind them and the tier names imply no confidence
interval. The SUPPORTED branch has never fired on real data, so the path that would escalate
a genuine notifiable positive is untested. Treat every verdict as triage, not determination,
and check positives by hand.
"""
from __future__ import annotations
import argparse, csv, json, re, sys
from pathlib import Path
from collections import defaultdict

# ------------------------------------------------- evidence thresholds (UNCALIBRATED)
# Calibrated against validations run on this cohort, where a refuted call and a confirmed
# call differed by roughly an order of magnitude on every axis:
#   refuted   S. pneumoniae : 3.9% of reads mapped,  7.6% genome breadth, dup 4.25
#   confirmed HCoV-OC43     : 100% of reads mapped, 99.9% genome breadth, identity 97.7%
DUP_SUSPECT      = 5.0    # KrakenUniq reads-per-unique-kmer above this looks like smearing
BREADTH_STRONG   = 50.0   # % of reference genome covered
BREADTH_WEAK     = 10.0
IDENTITY_REAL    = 92.0   # ONT vs a true species match; spurious alignment sits far below
IDENTITY_SPURIOUS= 85.0
MAPPED_FRAC_REAL = 50.0   # % of the reads assigned to a taxon that actually align to it


def norm(name: str) -> str:
    """Normalise a taxon label for comparison."""
    n = re.sub(r'\s+', ' ', (name or '')).strip()
    n = re.sub(r'\s*\(taxid \d+\)$', '', n)
    return n


def load_notifiable(path: Path):
    """Return (species -> [diseases]), (genus -> {species}), and genus->diseases."""
    species, genus_species, genus_disease = defaultdict(list), defaultdict(set), defaultdict(set)
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            g  = norm(row.get('Genus', ''))
            sp = norm(row.get('Species / taxon', ''))
            dis = norm(row.get('Disease / condition', ''))
            atype = norm(row.get('Agent type', ''))
            if g in ('', '—'):
                continue
            genus_disease[g].add(dis)
            if sp and sp != '—':
                species[sp].append((dis, atype, norm(row.get('Agent / strain / subgroup', ''))))
                genus_species[g].add(sp)
    return species, genus_species, genus_disease


def load_profile(path: Path):
    """Our per-sample species tables: <name>\t<taxid>\t<reads>\t<pct>."""
    out = []
    with path.open() as fh:
        rd = csv.reader(fh, delimiter='\t')
        hdr = next(rd, None)
        for r in rd:
            if len(r) < 3:
                continue
            try:
                reads = int(r[2]); pct = float(r[3]) if len(r) > 3 and r[3] not in ('', 'NA') else None
            except ValueError:
                continue
            out.append({'taxon': norm(r[0]), 'taxid': r[1], 'reads': reads, 'pct': pct})
    return out


def load_krakenuniq(path: Path):
    """taxName -> {reads, kmers, dup, cov} from a KrakenUniq report."""
    ev = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('%'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            try:
                ev[norm(f[8])] = {'reads': int(f[1]), 'taxreads': int(f[2]),
                                  'kmers': int(f[3]), 'dup': float(f[4]), 'cov': float(f[5])}
            except ValueError:
                continue
    return ev


def load_mapping(path: Path):
    """Optional mapping-breadth evidence: taxon\treference\treads_in\treads_mapped\tpct_mapped\tidentity\tbreadth."""
    ev = defaultdict(list)
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            key = norm(r.get('taxon') or r.get('sample') or '')
            ev[key].append(r)
    return ev


def assess(hit, kuniq, mapping, neg_taxa):
    """Return (confidence, verdict, reasons[]) - deliberately conservative."""
    reasons, score = [], 0
    name = hit['taxon']

    k = kuniq.get(name)
    if k:
        if k['dup'] > DUP_SUSPECT:
            score -= 2
            reasons.append(f"KrakenUniq duplication {k['dup']:.1f} (>{DUP_SUSPECT}) - few unique k-mers "
                           f"spread over many reads, the signature of cross-assignment")
        else:
            score += 1
            reasons.append(f"KrakenUniq duplication {k['dup']:.2f} with {k['kmers']:,} unique k-mers")
        if k['cov'] * 100 >= BREADTH_STRONG:
            score += 2; reasons.append(f"unique k-mers cover {k['cov']*100:.1f}% of the reference genome")
        elif k['cov'] * 100 < 0.01:
            # Near-zero genome coverage refutes a taxon even when duplication looks clean:
            # a handful of unique k-mers can sit at dup ~1 and still represent nothing.
            score -= 3
            reasons.append(f"unique k-mers cover only {k['cov']*100:.4f}% of the genome ({k['kmers']:,} k-mers) "
                           f"- far too little to support presence, regardless of duplication")
        elif k['cov'] * 100 < 1:
            score -= 1; reasons.append(f"unique k-mers cover only {k['cov']*100:.2f}% of the genome")
        # A large gap between the profile's read count and KrakenUniq's direct assignment is
        # itself evidence that the reads are being drawn in by shared sequence.
        if hit['reads'] > 50 and k['taxreads'] * 10 < hit['reads']:
            score -= 2
            reasons.append(f"profile reports {hit['reads']} reads but KrakenUniq assigns only "
                           f"{k['taxreads']} directly - most reads are not specific to this taxon")
    else:
        reasons.append("no KrakenUniq evidence available")

    # ---- discrimination against congeners -------------------------------------------
    # This is the single most informative test, and the one that settled every ambiguous
    # call in this cohort. A genuine population maps markedly better to its own reference
    # than to its relatives. Reads drawn in by conserved sequence map about equally well to
    # all of them, which is what we saw for S. pneumoniae, T. pallidum and N. meningitidis.
    refs = mapping.get(name, [])
    if len(refs) > 1:
        epithet = name.split()[-1].lower()
        own, others = [], []
        for m in refs:
            (own if epithet[:5] in m.get('reference', '').lower() else others).append(m)
        def br(m):
            try: return float(m.get('breadth_pct', 0) or 0)
            except ValueError: return 0.0
        if own and others:
            ob, xb = max(br(m) for m in own), max(br(m) for m in others)
            if xb > 0 and ob < xb * 1.2:
                score -= 4
                reasons.append(f"NO DISCRIMINATION: maps to its own reference ({ob:.1f}% breadth) no better "
                               f"than to congeners ({xb:.1f}%) - the signature of reads drawn in by "
                               f"sequence shared across the genus, not of this species being present")
            elif ob >= xb * 1.5:
                score += 3
                reasons.append(f"discriminates from congeners: {ob:.1f}% breadth on its own reference "
                               f"vs {xb:.1f}% on the nearest relative")

    best = None
    for m in mapping.get(name, []):
        try:
            br = float(m.get('breadth_pct', 0)); idt = float(m.get('mean_identity', 0) or 0)
        except ValueError:
            continue
        if best is None or br > best[0]:
            best = (br, idt, m)
    if best:
        br, idt, m = best
        pm = m.get('pct_mapped')
        if br >= BREADTH_STRONG: score += 3; reasons.append(f"maps across {br:.1f}% of the genome")
        elif br < BREADTH_WEAK:  score -= 2; reasons.append(f"maps across only {br:.1f}% of the genome")
        if idt:
            if idt >= IDENTITY_REAL:      score += 2; reasons.append(f"mean identity {idt:.1f}% - consistent with a true species match")
            elif idt < IDENTITY_SPURIOUS: score -= 3; reasons.append(f"mean identity {idt:.1f}% - too low for a real match; spurious alignment")
        if pm:
            try:
                pmf = float(pm)
                if pmf < MAPPED_FRAC_REAL:
                    score -= 3
                    reasons.append(f"only {pmf:.1f}% of the reads called this taxon actually align to its genome")
            except ValueError:
                pass
    else:
        reasons.append("no mapping-breadth evidence available")

    if name in neg_taxa:
        score -= 2; reasons.append(f"ALSO PRESENT IN THE NEGATIVE CONTROL ({neg_taxa[name]} reads) - possible contamination")

    if hit['reads'] < 10:
        score -= 1; reasons.append(f"only {hit['reads']} reads")

    have_validation = bool(k) or bool(best)
    if not have_validation:
        return 'UNVALIDATED', 'DO NOT REPORT - validation not yet run', reasons
    if score >= 4:   return 'STRONG',  'SUPPORTED - treat as possibly notifiable, escalate', reasons
    if score >= 1:   return 'MODERATE','UNRESOLVED - validate further before reporting', reasons
    if score >= -2:  return 'WEAK',    'UNRESOLVED - validate further before reporting', reasons
    return 'REFUTED', 'NOT SUPPORTED - do not report as present', reasons


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-n', '--notifiable', required=True, type=Path, help='notifiable organism TSV')
    ap.add_argument('-p', '--profile', required=True, type=Path, nargs='+', help='per-sample species table(s)')
    ap.add_argument('-k', '--krakenuniq', type=Path, help='KrakenUniq report for unique k-mer evidence')
    ap.add_argument('-m', '--mapping', type=Path, help='mapping-breadth TSV')
    ap.add_argument('--negative-control', type=Path, help='species table of a negative control')
    ap.add_argument('--min-reads', type=int, default=1)
    ap.add_argument('-o', '--output', type=Path, help='write TSV report here')
    ap.add_argument('--json', type=Path, help='also write JSON')
    args = ap.parse_args()

    species, genus_species, genus_disease = load_notifiable(args.notifiable)
    kuniq   = load_krakenuniq(args.krakenuniq) if args.krakenuniq else {}
    mapping = load_mapping(args.mapping) if args.mapping else {}
    neg_taxa = {}
    if args.negative_control and args.negative_control.exists():
        neg_taxa = {h['taxon']: h['reads'] for h in load_profile(args.negative_control) if h['reads'] > 0}

    findings = []
    for prof in args.profile:
        sample = prof.stem.replace('.species', '')
        for hit in load_profile(prof):
            if hit['reads'] < args.min_reads:
                continue
            name = hit['taxon']
            gen = name.split()[0] if name else ''
            if name in species:
                cat, diseases = 'SPECIES_NOTIFIABLE', sorted({d for d, _, _ in species[name]})
            elif gen in genus_species:
                unresolved = bool(re.search(r'\bsp\.?$|\bsp\.\s|\bspp\.?', name)) or name == gen
                cat = 'GENUS_UNRESOLVED' if unresolved else 'GENUS_CONGENER'
                diseases = sorted(genus_disease.get(gen, []))
            else:
                continue
            conf, verdict, reasons = assess(hit, kuniq, mapping, neg_taxa)
            if cat == 'GENUS_CONGENER':
                verdict = ('NOT NOTIFIABLE - different species in a notifiable genus; '
                           'record so the genus name is never reported bare')
            findings.append({'sample': sample, 'taxon': name, 'category': cat,
                             'reads': hit['reads'], 'pct': hit['pct'],
                             'notifiable_species_in_genus': sorted(genus_species.get(gen, [])),
                             'diseases': diseases, 'confidence': conf,
                             'verdict': verdict, 'evidence': reasons})

    order = {'SPECIES_NOTIFIABLE': 0, 'GENUS_UNRESOLVED': 1, 'GENUS_CONGENER': 2}
    findings.sort(key=lambda f: (order[f['category']], -f['reads']))

    # For every notifiable genus detected, state explicitly which notifiable species in that
    # genus were and were not found. An unstated absence is what caused the Treponema problem:
    # the genus was named, the notifiable species was never mentioned either way, and the
    # reader supplied the alarming interpretation. Absence must be asserted, not implied.
    detected = defaultdict(set)
    for prof in args.profile:
        sample = prof.stem.replace('.species', '')
        for hit in load_profile(prof):
            if hit['reads'] >= args.min_reads:
                detected[sample].add(hit['taxon'])
    explicit = []
    for f in findings:
        gen = f['taxon'].split()[0]
        if gen not in genus_species:
            continue
        for sp in sorted(genus_species[gen]):
            present = sp in detected[f['sample']]
            explicit.append({'sample': f['sample'], 'genus': gen, 'notifiable_species': sp,
                             'detected': present,
                             'statement': (f"{sp} DETECTED - see its own entry" if present else
                                           f"{sp} NOT detected in {f['sample']} (0 reads at the screening depth)")})
    seen_e = set(); explicit = [e for e in explicit
                                if not ((e['sample'], e['notifiable_species']) in seen_e
                                        or seen_e.add((e['sample'], e['notifiable_species'])))]

    if args.output:
        with args.output.open('w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['sample','taxon','category','reads','pct','confidence','verdict','diseases','evidence'])
            for f in findings:
                w.writerow([f['sample'],f['taxon'],f['category'],f['reads'],
                            '' if f['pct'] is None else f"{f['pct']:.4f}",
                            f['confidence'],f['verdict'],'; '.join(f['diseases']),' | '.join(f['evidence'])])
    if args.output:
        with Path(str(args.output).replace('.tsv', '') + '.explicit_negatives.tsv').open('w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['sample', 'genus', 'notifiable_species', 'detected', 'statement'])
            for e in explicit:
                w.writerow([e['sample'], e['genus'], e['notifiable_species'], e['detected'], e['statement']])
    if args.json:
        args.json.write_text(json.dumps({'findings': findings, 'explicit_statements': explicit}, indent=2))

    esc = [f for f in findings if f['category'] == 'SPECIES_NOTIFIABLE']
    amb = [f for f in findings if f['category'] == 'GENUS_UNRESOLVED']
    con = [f for f in findings if f['category'] == 'GENUS_CONGENER']
    print("=" * 78)
    print("NOTIFIABLE-ORGANISM SCREEN")
    print("=" * 78)
    print("  WARNING: confidence thresholds are NOT calibrated - set by eye from one cohort,")
    print("  no labelled truth set, no sensitivity/specificity estimate. The SUPPORTED branch")
    print("  has never fired on real data. Treat verdicts as TRIAGE, not determination, and")
    print("  verify any positive by hand. See pawsey_minion/AGENTS.md section 10.")
    print("-" * 78)
    print(f"  species-level matches to notifiable organisms : {len(esc)}")
    print(f"  unresolved calls in notifiable genera         : {len(amb)}   <- MUST resolve")
    print(f"  congeners (notifiable genus, other species)   : {len(con)}")
    blocking = [f for f in esc + amb if f['confidence'] in ('UNVALIDATED','MODERATE','WEAK','STRONG')
                and not f['verdict'].startswith('NOT SUPPORTED')]
    print(f"\n  ==> {len(blocking)} finding(s) must NOT be reported without the stated validation\n")
    for f in esc + amb:
        print(f"--- {f['sample']}  {f['taxon']}  [{f['category']}]")
        print(f"      reads={f['reads']}" + (f"  ({f['pct']:.3f}%)" if f['pct'] is not None else ''))
        if f['diseases']: print(f"      notifiable as: {', '.join(f['diseases'])}")
        print(f"      CONFIDENCE: {f['confidence']}   VERDICT: {f['verdict']}")
        for r in f['evidence']: print(f"        - {r}")
        print()
    print("=" * 78)
    print("EXPLICIT STATEMENTS FOR EVERY NOTIFIABLE SPECIES IN A DETECTED GENUS")
    print("(state these; never let an absence be inferred from silence)")
    print("=" * 78)
    for e in sorted(explicit, key=lambda x: (x['sample'], x['genus'], x['notifiable_species'])):
        mark = '!!' if e['detected'] else 'ok'
        print(f"  [{mark}] {e['sample']}: {e['statement']}")
    print()

    if con:
        print("Congeners recorded (report the SPECIES, never the bare genus):")
        seen = set()
        for f in con:
            key = (f['taxon'],)
            if key in seen: continue
            seen.add(key)
            print(f"    {f['taxon']:42} genus also contains: {', '.join(f['notifiable_species_in_genus'])}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
