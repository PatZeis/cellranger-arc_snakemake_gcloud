
configfile: 'config/config.yaml'

include: 'rules/common.smk'

wildcard_constraints:
    sample="[^\/]*"

rule results:
    input:
        cellranger = expand( '/data/cellranger_arc_run/{sample}.tar',
                     sample=samples_mk["id"].tolist(),
                     ),
            
    resources:
        machine_type = 'n1-standard-96', 
        disk_mb      = 10000000
             
### demultiplex bcl files for rna and atac layer and run cellranger 

rule create_simple_csv:
    input:
        'data/ref/mkfastq_samples.tsv' ### generate mkfastq_shample.tsv sheet using scripts/exectute_mksample_file.sh and store at google bucket in the folder bucket/google/ref
    output:
        'data/mkfastq/cellranger_arc_mkfastq_samplefile_{sample}.csv'
    params:
        sampleid = '{sample}'
                
    conda:
        'env/mkfastq_csv.yaml'
    
    resources:
        machine_type      = 'n1-standard-64',   
        disk_mb           = 4000000            
 
    shell:
        '''
        python src/create_mkfastq_simple_csv.py --samples {input.sample_file} --output_file {output} --wildcard {params.sampleid} 
        '''

rule demultiplex_bcl_file:
    input:
        cr = 'data/mkfastq/cellranger-arc-2.0.2.tar.gz', ### download cellranger archive from the 10x website and store at google bucket in bucket/data/mkfastq
        csv = 'data/mkfastq/cellranger_arc_mkfastq_samplefile_{sample}.csv',
        bcl = 'data/mkfastq/bcl_folder_rna.tar' ###  archive bcl folder: tar -cvf bcl_rna.tar bcl_rna
    output:
        mkfastq = 'data/mkfastq/{sample}.tar',  
        order_rna = 'data/mkfastq/{sample}_dump.tar'
    params:
        sampleid = '{sample}',  
        sampleid_tar = '{sample}.gz',
        bcl_folder = 'bcl_folder_rna',
        dump = '{sample}_dump',
        rna_fq_folder = expand( '{{sample}}/{fastq_path}',
                             fastq_path=config['fastq_path'],
                             ),
        reports = '{sample}/outs/fastq_path/Reports',
        stats = '{sample}/outs/fastq_path/Stats',
        ss = '{sample}/outs/input_samplesheet.csv'  
 
    conda:
        'env/mkfastq.yaml'
    resources:
        machine_type      = 'n1-highmem-96',
        disk_mb           = 20000000
    shell:
        '''
        tar -xzvf {input.cr}
        tar -xvf {input.bcl} 
        cellranger-arc-2.0.2/cellranger-arc mkfastq --id={params.sampleid} --run={params.bcl_folder} --csv={input.csv}
        rm -r {params.bcl_folder}
        mkdir temp_folder
        mv {params.rna_fq_folder} temp_folder/
        mv {params.reports} temp_folder/
        mv {params.stats} temp_folder/
        mv {params.ss} temp_folder/
        rm -r {params.sampleid}
        mv temp_folder {params.sampleid}
        tar -cvf {params.sampleid_tar} {params.sampleid}
        rm -r {params.sampleid}
        mv {params.sampleid_tar} {output.mkfastq}
        mkdir {params.dump}
        tar -cvf {output.order_rna} {params.dump}
        '''

rule create_simple_csv_atac:
    input:
        ordering = 'data/mkfastq/{sample}_dump.tar', 
        sample = 'data/ref/mkfastq_samples.tsv'
    output:
        'data/mkfastq/ATAC/cellranger_arc_mkfastq_samplefile_{sample}.csv'

    params:
        sampleid = '{sample}'

    conda:
        'env/mkfastq_csv.yaml'

    resources:
        machine_type      = 'n1-standard-64',   
        disk_mb           = 4000000             

    shell:
        '''
        python src/create_mkfastq_simple_csv.py --samples {input.sample} --output_file {output} --wildcard {params.sampleid} --atac go
        '''
        
rule demultiplex_bcl_file_atac:
    input:
        cr = 'data/mkfastq/cellranger-arc-2.0.2.tar.gz',
        csv = 'data/mkfastq/ATAC/cellranger_arc_mkfastq_samplefile_{sample}.csv',
        bcl = 'data/mkfastq/ATAC/bcl_folder_atac.tar'
    output:
        mkfastq = 'data/mkfastq/ATAC/{sample}.tar',  
        order_atac = 'data/mkfastq/ATAC/{sample}_dump.tar'
    params:
        sampleid = '{sample}',  
        sampleid_tar = '{sample}.gz',
        bcl_folder = 'bcl_folder_atac', 
        atac_fq_folder = expand( '{{sample}}/{fastq_path_atac}',
                             fastq_path_atac=config['fastq_path_atac'],
                             ),
        reports = '{sample}/outs/fastq_path/Reports',
        stats = '{sample}/outs/fastq_path/Stats',
        ss = '{sample}/outs/input_samplesheet.csv',
        dump = '{sample}_dump'
    conda:
        'env/mkfastq.yaml'
    resources:
        machine_type      = 'n1-highmem-96',
        disk_mb           = 20000000
    shell:
        '''
        tar -xzvf {input.cr}
        tar -xvf {input.bcl} 
        cellranger-arc-2.0.2/cellranger-arc mkfastq --id={params.sampleid} --run={params.bcl_folder} --csv={input.csv}
        rm -r {params.bcl_folder}
        mkdir temp_folder
        mv {params.atac_fq_folder} temp_folder/
        mv {params.reports} temp_folder/
        mv {params.stats} temp_folder/
        mv {params.ss} temp_folder/
        rm -r {params.sampleid}
        mv temp_folder {params.sampleid} 
        tar -cvf {params.sampleid_tar} {params.sampleid}
        rm -r {params.sampleid}
        mv {params.sampleid_tar} {output.mkfastq}
        mkdir {params.dump}
        tar -cvf {output.order_atac} {params.dump}
        '''


# generate library csv for cellranger arc

rule create_library_csv:
    input:
        mksample = 'data/ref/whitelist/mkfastq_samples.tsv',
        order1 = 'data/mkfastq/ATAC/{sample}_dump.tar'

    output:
        'data/mkfastq/ATAC/libraries_{sample}.csv'
    params:
        rna_path = expand( '/workdir/rna/{{sample}}/{flow_cell}',
                         flow_cell=config['flow_cell'],
                         ), 
        atac_path = expand( '/workdir/atac/{{sample}}/{flow_cell_atac}/{{sample}}_{fc}',
                         flow_cell_atac=config['flow_cell_atac'],
                         fc=config['fc'],
                         ),
        sampleid = '{sample}'
    conda:
        'env/mkfastq_csv.yaml'
    resources:
        machine_type      = 'n1-standard-64',   # vCPUs: 64, RAM: 60GB
        disk_mb           = 4000000             # 4 TB - needed for optimal read/write speed  
    shell:
        '''
        python src/create_cellranger_library.py --samples {input.mksample} --RNA_path {params.rna_path} --ATAC_path {params.atac_path}  --wildcard {params.sampleid} --output_file {output}
        '''

### run cellranger-arc count

rule cellranger_arc_count:
    input:
        cr = 'data/mkfastq/cellranger-arc-2.0.2.tar.gz',      
        ref_genome = 'data/ref/ref_genome_101_Macaca.gz', ### create reference genome according to 10x website, archive it and store on google bucket at bucket/data/ref
        rna = 'data/mkfastq/{sample}.tar',
        atac = 'data/mkfastq/ATAC/{sample}.tar',
        library = 'data/mkfastq/ATAC/libraries_{sample}.csv'
    output:
        cr = '/data/cellranger_arc_run/{sample}.tar',
        order_cr = '/data/cellranger_arc_run/{sample}_dump.tar'
    params:
        sampleid = '{sample}',
        dump = '{sample}_dump'
    conda:
        'env/mkfastq.yaml'
    resources:
        machine_type      = 'n1-highmem-96',
        disk_mb           = 20000000
    shell:
        '''
        tar -xzvf {input.cr}
        mkdir rna
        mkdir atac
        mkdir reference
        tar -xvf {input.rna} -C ./rna
        tar -xvf {input.atac} -C ./atac
        tar -xzvf {input.ref_genome} -C ./reference --strip-components=1
        cellranger-arc-2.0.2/cellranger-arc count --id={params.sampleid} --reference=./reference --libraries={input.library} 
        tar -cvf {output.cr} {params.sampleid}
        mkdir {params.dump}
        tar -cvf {output.order_cr} {params.dump}
        '''      

        
        
