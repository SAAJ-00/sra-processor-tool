V0.2

Instrucciones de instalación:
  1. Clonar repositorio o descargar zip
     
  1.1 En caso de descargar zip descomprimir
  
  2. Crear ambiente conda con el archivo environment.yml
  
    conda env create -f environment.yml
  
  3. Activar ambiente
  
    conda activate sra_processor
  
  4. Instalar programa en modo editable desde carpeta principal
  
    pip install -e .
  
  5. Probar programa
  
    sra-processor -h
