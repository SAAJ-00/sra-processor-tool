V0.2

Instrucciones de instalación:
  1. Clonar repositorio o descargar zip
    1.1 En caso de descargar zip descomprimir
  3. Crear ambiente conda con el archivo ´´´environment.yml´´´
     $ conda env create -f environment.yml
  4. Activar ambiente
     $ conda activate sra_processor
  5. Instalar programa en modo editable desde carpeta principal
     $ pip install -e .
  6. Probar programa
     $ sra-processor -h
