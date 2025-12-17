# SRA Processor Tool

V0.3 - Herramienta modular para procesamiento de datos desde SRA IDs (Illumina/Nanopore/PacBio)

## Instalación

### 1. Clonar repositorio
```bash
git clone <repository-url>
cd sra-processor-tool
```

### 2. Crear ambiente conda
```bash
conda env create -f environment.yml
```

### 3. Activar ambiente
```bash
conda activate sra_processor
```

### 4. Instalar programa
```bash
pip install -e .
```

### 5. Verificar instalación
```bash
sra-processor -h
```

## Uso

La herramienta ofrece tres subcomandos principales para mayor flexibilidad:

### Subcomandos Disponibles

#### 1. `download` - Solo descarga y conversión
Descarga archivos SRA y los convierte a FASTQ sin hacer trimming.

```bash
sra-processor download <srr_list.txt> [opciones]
```

**Opciones:**
- `-o, --output`: Directorio de salida (default: `.`)
- `-t, --threads`: Número de hilos (default: 8)
- `--max-size`: Tamaño máximo por SRR (default: `30G`)
- `--keep-sra`: Mantener archivos .sra después de conversión
- `--keep-temp`: Mantener archivos temporales

**Ejemplo:**
```bash
# Descargar y convertir, manteniendo archivos SRA
sra-processor download srr_list.txt -o data/ --keep-sra

# Descargar con 16 hilos
sra-processor download srr_list.txt -o data/ -t 16 --max-size 100G
```

#### 2. `trim` - Solo trimming
Procesa archivos FASTQ existentes (de SRRs descargados o archivos externos).

```bash
sra-processor trim <input_file.txt> [opciones]
```

**Opciones:**
- `-o, --output`: Directorio de salida (default: `.`)
- `-t, --threads`: Número de hilos (default: 8)
- `--input-type`: Tipo de input: `auto`, `srr`, o `fastq` (default: `auto`)
- `--force`: Sobrescribir archivos existentes
- `--quality-phred`: Calidad Phred mínima (default: 30)
- `--min-length`: Longitud mínima de reads (default: 50 bp)
- Opciones adicionales para lecturas largas (ver `--help`)

**Ejemplos:**
```bash
# Trimming de SRRs ya descargados
sra-processor trim srr_list.txt -o data/ --quality-phred 30

# Trimming de archivos FASTQ externos
sra-processor trim fastq_files.txt -o results/ --input-type fastq

# Forzar re-procesamiento con calidad más estricta
sra-processor trim srr_list.txt -o data/ --force --quality-phred 35
```

#### 3. `full` - Pipeline completo
Ejecuta el pipeline completo: descarga, conversión y trimming.

```bash
sra-processor full <srr_list.txt> [opciones]
```

**Opciones:** Combina todas las opciones de `download` y `trim`.

**Ejemplo:**
```bash
# Pipeline completo con configuración estándar
sra-processor full srr_list.txt -o data/ -t 16

# Pipeline con parámetros personalizados
sra-processor full srr_list.txt -o data/ -t 16 \
  --max-size 50G \
  --quality-phred 35 \
  --min-length 75
```

## Casos de Uso Comunes

### 1. Descarga en Batch
Descargar múltiples SRRs sin procesar:
```bash
sra-processor download srr_list.txt -o raw_data/ --keep-sra -t 16
```

### 2. Procesar FASTQs Externos
Hacer trimming de archivos FASTQ que ya tienes:
```bash
# Crear archivo con rutas a FASTQs
echo "/path/to/sample1_1.fastq" > fastq_list.txt
echo "/path/to/sample1_2.fastq" >> fastq_list.txt
echo "/path/to/sample2_1.fastq" >> fastq_list.txt

# Procesar
sra-processor trim fastq_list.txt -o trimmed/ --input-type fastq
```

### 3. Reanudar Pipeline Interrumpido
Si el pipeline se interrumpe, simplemente vuelve a ejecutar el mismo comando:
```bash
# La herramienta detecta automáticamente el estado y reanuda
sra-processor full srr_list.txt -o data/ -t 16
```

### 4. Re-procesar con Parámetros Diferentes
Para re-hacer trimming con parámetros diferentes:
```bash
sra-processor trim srr_list.txt -o data/ --force \
  --quality-phred 35 \
  --min-length 100
```

## Estados del Sistema

La herramienta detecta automáticamente el estado de cada SRR:

| Estado | Descripción | Acción al reanudar |
|--------|-------------|-------------------|
| `new` | No iniciado | Inicia desde cero |
| `sra_downloaded` | SRA descargado | Convierte a FASTQ |
| `fastq_ready` | FASTQ generados (paired-end) | Hace trimming |
| `single_end` | FASTQ single-end generado | Hace trimming |
| `complete` | Procesamiento completo | Omite (usa `--force` para re-procesar) |
| `unknown` | Estado no reconocido | Limpia y reinicia |

## Formato de Archivos de Entrada

### Lista de SRRs (`srr_list.txt`)
```
SRR1234567
SRR1234568
SRR1234569
```

### Lista de FASTQs (`fastq_list.txt`)
```
/path/to/sample1_1.fastq
/path/to/sample1_2.fastq
/path/to/sample2.fastq
/data/sample3_1.fastq.gz
/data/sample3_2.fastq.gz
```

## Estructura de Directorios de Salida

```
output_dir/
├── SRR1234567/
│   ├── SRR1234567_1.fastq          # FASTQ raw (paired-end)
│   ├── SRR1234567_2.fastq
│   ├── SRR1234567_1_trimmed.fastq.gz  # FASTQ procesado
│   ├── SRR1234567_2_trimmed.fastq.gz
│   └── report_SRR1234567.html      # Reporte de calidad
├── SRR1234568/
│   ├── SRR1234568.fastq            # FASTQ raw (single-end)
│   ├── SRR1234568_trimmed.fastq.gz # FASTQ procesado
│   └── report_SRR1234568.html
└── logs/
    └── sra_processor.log           # Log completo
```

## Tipos de Datos Soportados

- **Illumina (paired-end)**: Detectado automáticamente
- **Illumina (single-end)**: Detectado automáticamente
- **Nanopore/PacBio (long reads)**: Detectado por longitud promedio (>500bp) o usar `--force-long-reads`

## Opciones Avanzadas

### Para lecturas cortas (Illumina)
```bash
sra-processor trim srr_list.txt -o data/ \
  --quality-phred 30 \
  --min-length 50 \
  --cut-window-size 4 \
  --cut-mean-quality 25
```

### Para lecturas largas (Nanopore/PacBio)
```bash
sra-processor trim srr_list.txt -o data/ \
  --force-long-reads \
  --long-min-quality 10 \
  --long-min-length 1000
```

### Deshabilitar filtros específicos
```bash
sra-processor trim srr_list.txt -o data/ \
  --disable-adapter-trimming \
  --disable-quality-filtering
```

## Troubleshooting

### Error: "Herramientas SRA no encontradas"
Asegúrate de tener instalado SRA Toolkit:
```bash
conda install -c bioconda sra-tools
```

### Error: "fastp no encontrado"
Instala fastp:
```bash
conda install -c bioconda fastp
```

### Para lecturas largas, instala fastplong:
```bash
# Instalación según documentación de fastplong
```

## Changelog

### V0.3
- ✨ Nueva arquitectura modular con subcomandos
- ✨ Subcomando `download`: solo descarga y conversión
- ✨ Subcomando `trim`: solo trimming (soporta SRRs y FASTQs externos)
- ✨ Subcomando `full`: pipeline completo (comportamiento anterior)
- ✨ Opción `--keep-sra` para mantener archivos SRA
- ✨ Opción `--force` para sobrescribir archivos existentes
- ✨ Auto-detección de tipo de input (SRR vs FASTQ)
- 🔧 Métodos públicos en SRADownloader para mayor flexibilidad
- 🔧 Soporte mejorado para archivos FASTQ externos en trimmer
- 📚 Documentación ampliada con ejemplos de uso

### V0.2
- Pipeline básico funcional

## Licencia

[Especificar licencia]

## Contribuciones

[Especificar cómo contribuir]

## Contacto

[Información de contacto]
