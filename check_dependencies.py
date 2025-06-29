# check_dependencies.py
import shutil
import sys
from pathlib import Path

REQUIRED_TOOLS = {
    'prefetch': 'Herramienta de SRA Toolkit para descarga',
    'fasterq-dump': 'Herramienta de SRA Toolkit para conversión',
    'fastp': 'Herramienta para procesamiento de FASTQ',
    'fastplong': 'Herramienta para procesamiento de lecturas largas',
    'pigz': 'Compresión paralela de archivos'
}

def check_tools():
    missing = []
    for tool, description in REQUIRED_TOOLS.items():
        if not shutil.which(tool):
            missing.append((tool, description))
    
    if missing:
        print("ERROR: Faltan las siguientes herramientas requeridas:", file=sys.stderr)
        for tool, desc in missing:
            print(f"- {tool}: {desc}", file=sys.stderr)
        
        print("\nPuedes instalarlas con:", file=sys.stderr)
        print("conda install -c bioconda sra-tools fastp fastplong pigz", file=sys.stderr)
        sys.exit(1)
    else:
        print("Todas las dependencias están instaladas correctamente")

if __name__ == '__main__':
    check_tools()