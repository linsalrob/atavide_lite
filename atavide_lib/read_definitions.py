"""
Read the definitions file and return a dictionary of the definitions
"""

def read_definitions(defintions_file = "DEFINITIIONS.sh"):
    """
    Read the definitions file and return a dictionary of the definitions
    """

    definitions = {}
    with open(defintions_file, 'r') as f:
        for l in f:
            if l.startswith("export"):
                p = l.replace("export ", "").strip().split("=")
                definitions[p[0]] = p[1].replace('"', '')
    return definitions