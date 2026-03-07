import subprocess
import tempfile
import os
import io
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def run_al2co(msa: MultipleSeqAlignment, executable: str = './al2co', **kwargs):
    """
    Run al2co on a BioPython MultipleSeqAlignment object and return the conservation scores.
    """
    sio = io.StringIO()
    AlignIO.write(msa, sio, "clustal")
    content = sio.getvalue()
    if content.startswith("CLUSTAL X"):
        content = "CLUSTAL W" + content[9:]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as f_in:
        f_in.write(content)
        temp_in = f_in.name

    temp_out = temp_in + ".out"

    args = [executable, '-i', temp_in, '-o', temp_out]
    for k, v in kwargs.items():
        if isinstance(v, bool):
            if v:
                args.extend(['-' + k, 'T'])
            else:
                args.extend(['-' + k, 'F'])
        else:
            args.extend(['-' + k, str(v)])

    try:
        result = subprocess.run(args, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        os.remove(temp_in)
        if os.path.exists(temp_out):
            os.remove(temp_out)
        raise RuntimeError(f"al2co failed: {e.stderr}\n{e.stdout}")

    scores = {}
    with open(temp_out, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*') or line.startswith('M:') or line.startswith('al2co') or line.startswith('Input') or line.startswith('Output') or line.startswith('Conservation') or line.startswith('Gap') or line.startswith('Window') or line.startswith('Weighting'):
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                pos = int(parts[0])
                score = float(parts[2])
                scores[pos] = score

    os.remove(temp_in)
    os.remove(temp_out)
    return scores
