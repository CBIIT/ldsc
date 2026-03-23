import os
import glob
import subprocess
import shutil
import pandas as pd
import gzip
from ldscore import parse as ps
from ldscore import ldscore as ld


def validSumstats(sumstats_file):
    """
    Validate summary statistics file format before processing.
    
    Args:
        sumstats_file: Path to the summary statistics file
        
    Returns:
        dict: Validation result with keys:
            - 'valid': bool indicating if file is valid
            - 'errors': list of error messages
            - 'warnings': list of warning messages
            - 'columns': list of detected columns
            - 'mapped_columns': dict of column mappings
    """
    result = {
        'valid': True,
        'errors': [],
        'warnings': [],
        'columns': [],
        'mapped_columns': {}
    }
    
    # Default column name mappings (from munge_sumstats.py)
    default_cnames = {
        'SNP': 'SNP', 'MARKERNAME': 'SNP', 'SNPID': 'SNP', 'RS': 'SNP', 'RSID': 'SNP',
        'RS_NUMBER': 'SNP', 'RS_NUMBERS': 'SNP',
        'P': 'P', 'PVALUE': 'P', 'P_VALUE': 'P', 'PVAL': 'P', 'P_VAL': 'P', 'GC_PVALUE': 'P',
        'A1': 'A1', 'ALLELE1': 'A1', 'ALLELE_1': 'A1', 'EFFECT_ALLELE': 'A1',
        'REFERENCE_ALLELE': 'A1', 'INC_ALLELE': 'A1', 'EA': 'A1',
        'A2': 'A2', 'ALLELE2': 'A2', 'ALLELE_2': 'A2', 'OTHER_ALLELE': 'A2',
        'NON_EFFECT_ALLELE': 'A2', 'DEC_ALLELE': 'A2', 'NEA': 'A2',
        'N': 'N', 'NCASE': 'N_CAS', 'CASES_N': 'N_CAS', 'N_CASE': 'N_CAS',
        'N_CASES': 'N_CAS', 'N_CONTROLS': 'N_CON', 'N_CAS': 'N_CAS', 'N_CON': 'N_CON',
        'NCONTROL': 'N_CON', 'CONTROLS_N': 'N_CON', 'N_CONTROL': 'N_CON',
        'ZSCORE': 'Z', 'Z-SCORE': 'Z', 'GC_ZSCORE': 'Z', 'Z': 'Z',
        'OR': 'OR', 'B': 'BETA', 'BETA': 'BETA', 'LOG_ODDS': 'LOG_ODDS',
        'EFFECTS': 'BETA', 'EFFECT': 'BETA',
        'INFO': 'INFO',
        'EAF': 'FRQ', 'FRQ': 'FRQ', 'MAF': 'FRQ', 'FRQ_U': 'FRQ', 'F_U': 'FRQ'
    }
    
    try:
        # Check if file exists
        if not os.path.exists(sumstats_file):
            result['valid'] = False
            result['errors'].append(f"File not found: {sumstats_file}")
            return result
        
        # Read header
        try:
            if sumstats_file.endswith('.gz'):
                with gzip.open(sumstats_file, 'rt') as f:
                    header_line = f.readline()
            else:
                with open(sumstats_file, 'r') as f:
                    header_line = f.readline()
            
            columns = header_line.strip().split()
            result['columns'] = columns
            
            if len(columns) == 0:
                result['valid'] = False
                result['errors'].append("File has no columns in header")
                return result
                
        except Exception as e:
            result['valid'] = False
            result['errors'].append(f"Error reading file header: {str(e)}")
            return result
        
        # Map column names (case-insensitive)
        def clean_header(s):
            return s.upper().replace('_', '').replace('.', '')
        
        for col in columns:
            clean_col = clean_header(col)
            if clean_col in default_cnames:
                result['mapped_columns'][col] = default_cnames[clean_col]
        
        # Check for required columns
        mapped_values = set(result['mapped_columns'].values())
        
        # Must have SNP
        if 'SNP' not in mapped_values:
            result['valid'] = False
            result['errors'].append("Missing required column: SNP (variant identifier)")
        
        # Must have P
        if 'P' not in mapped_values:
            result['valid'] = False
            result['errors'].append("Missing required column: P (p-value)")
        
        # Must have at least one signed statistic (Z, OR, BETA, or LOG_ODDS)
        signed_stats = {'Z', 'OR', 'BETA', 'LOG_ODDS'}
        if not any(stat in mapped_values for stat in signed_stats):
            result['valid'] = False
            result['errors'].append("Missing signed summary statistic column (need one of: Z, OR, BETA, LOG_ODDS)")
        
        # Check for sample size column - now treated as error
        if 'N' not in mapped_values and 'N_CAS' not in mapped_values and 'N_CON' not in mapped_values:
            result['valid'] = False
            result['errors'].append("No sample size column found (N, N_CAS, N_CON). You must provide --N or --N-cas/--N-con")
        
        # Check for allele columns (recommended but not required)
        if 'A1' not in mapped_values or 'A2' not in mapped_values:
            result['warnings'].append("Missing allele columns (A1, A2). This is OK for h2 estimation but required for genetic correlation")
        
        # Try to read a few rows to validate data types
        try:
            if sumstats_file.endswith('.gz'):
                df_sample = pd.read_csv(sumstats_file, sep=r'\s+', nrows=10, compression='gzip')
            else:
                df_sample = pd.read_csv(sumstats_file, sep=r'\s+', nrows=10)
            
            # Check if numeric columns are numeric - now treated as error
            for col in df_sample.columns:
                if col in result['mapped_columns']:
                    mapped = result['mapped_columns'][col]
                    if mapped in ['P', 'N', 'N_CAS', 'N_CON', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ']:
                        if not pd.api.types.is_numeric_dtype(df_sample[col]):
                            result['valid'] = False
                            result['errors'].append(f"Column '{col}' (mapped to {mapped}) must be numeric")
                            
        except Exception as e:
            result['warnings'].append(f"Could not validate data types: {str(e)}")
        
    except Exception as e:
        result['valid'] = False
        result['errors'].append(f"Unexpected error during validation: {str(e)}")
    
    return result


def validBfile(bfile_prefix):
    """
    Validate PLINK bfile format (.bed/.bim/.fam) before processing.
    
    Args:
        bfile_prefix: Path prefix for PLINK files (without extension)
        
    Returns:
        dict: Validation result with keys:
            - 'valid': bool indicating if files are valid
            - 'errors': list of error messages
            - 'warnings': list of warning messages
            - 'n_snps': number of SNPs in .bim file (if valid)
            - 'n_samples': number of samples in .fam file (if valid)
    """
    result = {
        'valid': True,
        'errors': [],
        'warnings': [],
        'n_snps': None,
        'n_samples': None
    }
    
    try:
        # Define file paths
        bim_file = bfile_prefix + '.bim'
        fam_file = bfile_prefix + '.fam'
        bed_file = bfile_prefix + '.bed'
        
        # Check if all three files exist
        missing_files = []
        if not os.path.exists(bim_file):
            missing_files.append(f"{bim_file} (.bim file)")
        if not os.path.exists(fam_file):
            missing_files.append(f"{fam_file} (.fam file)")
        if not os.path.exists(bed_file):
            missing_files.append(f"{bed_file} (.bed file)")
        
        if missing_files:
            result['valid'] = False
            result['errors'].append(f"Missing required files: {', '.join(missing_files)}")
            return result
        
        # Validate .bim file
        try:
            # Check if first row might be a header by comparing data patterns with second row
            with open(bim_file, 'r') as f:
                lines = [f.readline().strip() for _ in range(2)]
                if len(lines) >= 2 and lines[0] and lines[1]:
                    row1_parts = lines[0].split()
                    row2_parts = lines[1].split()
                    
                    # .bim files have 6 columns: CHR SNP CM BP A1 A2
                    if len(row1_parts) >= 1 and len(row2_parts) >= 1:
                        # Compare column 1 (CHR) - should be numeric (1-22) or 'X', 'Y', 'MT'
                        chr1 = row1_parts[0]
                        chr2 = row2_parts[0]
                        
                        # Heuristic: Row 1 is likely a header if:
                        # - CHR column is alphabetic (like "CHR", "CHROM", "chromosome")
                        # - Row 2 CHR is numeric or valid chromosome value
                        is_row1_alpha = chr1.isalpha() and chr1.upper() not in ['X', 'Y', 'MT', 'XY', 'PAR1', 'PAR2']
                        is_row2_valid_chr = chr2.isdigit() or chr2.upper() in ['X', 'Y', 'MT', 'XY', 'PAR1', 'PAR2']
                        
                        if is_row1_alpha and is_row2_valid_chr:
                            result['warnings'].append(
                                f"BIM file appears to have a header. Column 1 (CHR) has '{chr1}' (looks like header) "
                                f"while second row has '{chr2}' (looks like data). "
                                f"PLINK .bim files should not have headers. "
                                f"The header will be treated as a variant, which may cause incorrect results."
                            )
            
            bim = ps.PlinkBIMFile(bim_file)
            result['n_snps'] = bim.n
            
            # Validate that numeric columns contain numeric values
            # CHR (column 0), CM (column 2), BP (column 3) should be numeric
            if 'CM' in bim.df.columns:
                if not pd.api.types.is_numeric_dtype(bim.df['CM']):
                    result['valid'] = False
                    result['errors'].append(
                        f"Column 'CM' (centimorgan) in .bim file must contain numeric values. "
                        f"Found non-numeric data which will cause errors during LD score calculation."
                    )
            
            if 'BP' in bim.df.columns:
                if not pd.api.types.is_numeric_dtype(bim.df['BP']):
                    result['valid'] = False
                    result['errors'].append(
                        f"Column 'BP' (base pair position) in .bim file must contain numeric values. "
                        f"Found non-numeric data."
                    )
        except ValueError as e:
            result['valid'] = False
            result['errors'].append(f"Invalid .bim file: {str(e)}")
        except Exception as e:
            result['valid'] = False
            result['errors'].append(f"Error reading .bim file: {str(e)}")
        
        # Validate .fam file
        try:
            # Check if first row might be a header by comparing data patterns with second row
            with open(fam_file, 'r') as f:
                lines = [f.readline().strip() for _ in range(2)]
                if len(lines) >= 2 and lines[0] and lines[1]:
                    row1_parts = lines[0].split()
                    row2_parts = lines[1].split()
                    
                    # Only check if both rows have at least 2 columns
                    if len(row1_parts) >= 2 and len(row2_parts) >= 2:
                        # Compare column 2 (IID) - LDSC only uses this column
                        iid1 = row1_parts[1]
                        iid2 = row2_parts[1]
                        
                        # Heuristic: Row 1 is likely a header if:
                        # - It's purely alphabetic (like "IID", "ID", "SAMPLE")
                        # - Row 2 contains digits or special chars (like "id1", "sample_001", "HG00096")
                        is_row1_alpha_only = iid1.isalpha()
                        is_row2_has_digits = any(c.isdigit() for c in iid2) or '_' in iid2 or '-' in iid2
                        
                        if is_row1_alpha_only and is_row2_has_digits:
                            result['warnings'].append(
                                f"First row appears to be a header. Column 2 (IID) has '{iid1}' (alphabetic only) "
                                f"while second row has '{iid2}' (contains digits/special chars). "
                                f"PLINK .fam files should not have headers. "
                                f"The header will be treated as a sample, which may cause incorrect results."
                            )
            
            # LDSC only reads column 2 (IID) with usecols=[1]
            fam = ps.PlinkFAMFile(fam_file)
            result['n_samples'] = fam.n
        except ValueError as e:
            result['valid'] = False
            result['errors'].append(f"Invalid .fam file: {str(e)}")
        except Exception as e:
            result['valid'] = False
            result['errors'].append(f"Error reading .fam file: {str(e)}")
        
        # Validate .bed file (only if .bim and .fam are valid)
        if result['n_snps'] is not None and result['n_samples'] is not None:
            try:
                bed = ld.PlinkBEDFile(bed_file, result['n_samples'], bim)
                # If we get here, the bed file is valid
                if bed.m == 0:
                    result['warnings'].append("No polymorphic SNPs remain after MAF filtering")
                elif bed.m < result['n_snps']:
                    result['warnings'].append(f"Only {bed.m} of {result['n_snps']} SNPs are polymorphic (monomorphic SNPs removed)")
            except IOError as e:
                result['valid'] = False
                if "Magic number" in str(e):
                    result['errors'].append(f"Invalid .bed file: Magic number not recognized. File may be corrupted or not a valid PLINK .bed file")
                elif "SNP-major mode" in str(e):
                    result['errors'].append(f"Invalid .bed file: Must be in SNP-major mode (default PLINK format)")
                elif "bits, expected" in str(e):
                    # Extract the actual and expected bit counts from the error message
                    result['errors'].append(
                        f"File size mismatch: .bed file size doesn't match .bim/.fam files. "
                        f"This usually indicates inconsistent files or corrupted .fam file. "
                        f"Expected {result['n_samples']} samples from .fam file, but .bed file size suggests a different number. {str(e)}"
                    )
                else:
                    result['errors'].append(f"Invalid .bed file: {str(e)}")
            except ValueError as e:
                result['valid'] = False
                result['errors'].append(f"Invalid .bed file: {str(e)}")
            except Exception as e:
                result['valid'] = False
                result['errors'].append(f"Error reading .bed file: {str(e)}")
        
    except Exception as e:
        result['valid'] = False
        result['errors'].append(f"Unexpected error during validation: {str(e)}")
    
    return result


def run_ldsc_command(pop, genome_build, filename,ldwindow,windUnit,isExample,reference):
    fileDir = f"/data/tmp/uploads/"
    if isinstance(isExample, str):
            isExample = isExample.lower() == 'true'
  
    
    #if isExample:
    #    fileDir =  "/data/ldscore"
    ldwindow_value = 1  # Example value, replace with actual value
    # Check if ldwindow is an integer greater than 0, if not set it to 1
    try:
        ldwindow_value = int(ldwindow)
        if ldwindow_value <= 0:
            ldwindow_value = 1
    except ValueError:
        ldwindow_value = 1

    windFlag = '--ld-wind-cm'
    if windUnit == 'cm':
        windFlag = "--ld-wind-cm"
    elif windUnit == 'kb':
        windFlag = "--ld-wind-kb"

    if filename:
        #file_parts = filename.split('.')
        file_chromo = filename
        # for part in file_parts:
        #     if part.isdigit() and 1 <= int(part) <= 22:
        #         file_chromo = part
        #         break
    
    # if file_chromo:
    #     # Find the file in the directory
    #     pattern = os.path.join(fileDir, f"{filename}.*")
    #     for file_path in glob.glob(pattern):
    #         extension = file_path.split('.')[-1]
    #         new_filename = f"{file_chromo}.{extension}"
    #         new_file_path = os.path.join(fileDir, new_filename)
    #         #os.rename(file_path, new_file_path)
    #         shutil.copy(file_path, new_file_path)  # Copy the file instead of renaming it
    
    file_chr=file_chromo
    if isExample:
        file_chr =  "/data/ldscore/"+file_chromo 
        fileDir = f"/data/tmp/uploads/{reference}/"  
    else:
        fileDir = f"/data/tmp/uploads/{reference}/"   
    try:
        # Run the command
        # 'cd 1kg_eur && python ../ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22'
        parent_dir = '/usr/local/bin/'
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print(fileDir)
        command = f"cd {fileDir} && python3 {ldsc_script_path} --bfile {file_chr} --l2 {windFlag} {ldwindow_value}  --out {file_chromo}"
        result = subprocess.run(
            ['bash', '-c', command],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"An error occurred: {e.stderr}"
    

#def main():
    # Example parameters for testing
#    pop = 'example_pop'
#    genome_build = 'example_genome_build'
#    filename = 'example_filename.22.txt'
    
    # Call the function and print the result
#    result = run_ldsc_command(pop, genome_build, filename)
#    print(result)

#if __name__ == "__main__":
#    main()

def run_herit_command(sumstats_file, fileDir, ld_scores_dir,isExample):
    fallExampleDir = f"/data/ldscore"
    w_hm3_snplist = f"/data/ldscore/w_hm3.snplist"
    if isinstance(isExample, str):
        isExample = isExample.lower() == 'true'
    errormsg = ""
    try:
        parent_dir = '/usr/local/bin/'
        munge_sumstat_script_path = os.path.join(parent_dir, 'munge_sumstats.py')
        # Generate the output filename based on the input summary statistics file
        base_name = os.path.splitext(os.path.basename(sumstats_file))[0]
        out_file = f"{base_name}.sumstats.gz"
       
                # If isExample is True, use the fallbackDir
        if isExample:
            sumstats_path = os.path.join(fallExampleDir, sumstats_file)
        else:
            sumstats_path = sumstats_file
        print("First command ################:",sumstats_path, isExample)
                # Ensure ld_scores_dir is in lowercase
        ld_scores_dir = ld_scores_dir.lower()
        ld_scores_dir = fallExampleDir+"/"+ld_scores_dir
        # Ensure ld_scores_dir has a trailing slash
        if not ld_scores_dir.endswith('/'):
            ld_scores_dir += '/'
        # First command
        command1 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path} --merge-alleles {w_hm3_snplist}  --out {base_name}"
        
        #command1 = f"python ../munge_sumstats.py --sumstats {sumstats_path} --merge-alleles ../testData/w_hm3.snplist  --out {base_name}"
      
        result1 = subprocess.run( ['bash', '-c', command1], check=True, capture_output=True, text=True)
       
        # command1 = [
        #     'python', '../munge_sumstats.py',
        #     '--sumstats', sumstats_file,
        #     '--merge-alleles', '../testData/w_hm3.snplist',
        #     '--a1', 'ALT',
        #     '--a2', 'REF',
        #     '--chunksize', '500000',
        #     '--out', base_name
        # ]
        #result1 = subprocess.run( command1, check=True, capture_output=True, text=True)
       
        print("First command output:", result1.stdout)
        #print("First command error (if any):", result1.stderr)

        # Second command
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print("Second command ################:",ld_scores_dir)
        command2 = f"cd {fileDir} && python3 {ldsc_script_path} --h2 {out_file} --ref-ld-chr {ld_scores_dir} --w-ld-chr {ld_scores_dir} --out {base_name}"
       
        #result2 = subprocess.run( ['bash', '-c', command2], check=True, capture_output=True, text=True)
        try:
            result2 = subprocess.run(['bash', '-c', command2], check=True, capture_output=True, text=True)
            print("Second command output:", result2.stdout, result2.stderr)
            separator = "\n--------\n"
            return result1.stdout + separator + result2.stdout
        except subprocess.CalledProcessError as e:
            error_details = (e.stderr or e.stdout or e.output or str(e)).strip()
            print(f"An error occurred while running the second command: {command2}")
            print(f"Second command exit code: {e.returncode}")
            print(f"Command stdout: {e.stdout}")
            print(f"Command stderr: {e.stderr}")
            separator = "\n--------\n"
            return f"{result1.stdout}{separator}Second command failed (exit {e.returncode}).\n{error_details}"

        # command2 = [
        #     'python', '../ldsc.py',
        #     '--h2', "test2.sumstats.gz",
        #     '--ref-ld-chr', ld_scores_dir,
        #     '--w-ld-chr', ld_scores_dir,
        #     '--out', base_name
        # ]
        # result2 = subprocess.run( command2, check=True, capture_output=True, text=True)
       
        #print("Second command output:", result2.stdout)
        #separator = "\n--------\n"
        #return result1.stdout + separator + errormsg
        #print("Second command error (if any):", result2.stderr)

    except subprocess.CalledProcessError as e:
        error_details = (e.stderr or e.stdout or e.output or str(e)).strip()
        print(f"An error occurred while running the command: {e}")
        print(f"Command output: {e.output}")
        print(f"Command stderr: {e.stderr}")
        return f"An error occurred while running the command (exit {e.returncode}): {error_details}"

def run_correlation_command(sumstats_file, sumstats_file2, fileDir, ld_scores_dir,isExample):
    fallExampleDir = f"/data/ldscore"
    w_hm3_snplist = f"/data/ldscore/w_hm3.snplist"
    if isinstance(isExample, str):
        isExample = isExample.lower() == 'true'
    errormsg = ""
    try:
        parent_dir = '/usr/local/bin/'
        munge_sumstat_script_path = os.path.join(parent_dir, 'munge_sumstats.py')
        # Generate the output filename based on the input summary statistics file
        base_name = os.path.splitext(os.path.basename(sumstats_file))[0]
        out_file = f"{base_name}.sumstats.gz"
        base_name2 = os.path.splitext(os.path.basename(sumstats_file2))[0]
        out_file2 = f"{base_name2}.sumstats.gz"
       
                # If isExample is True, use the fallbackDir
        if isExample:
            sumstats_path = os.path.join(fallExampleDir, sumstats_file)
            sumstats_path2 = os.path.join(fallExampleDir, sumstats_file2)
        else:
            sumstats_path = sumstats_file
            sumstats_path2 = sumstats_file2
        print("First command ################:",sumstats_path, isExample)
                # Ensure ld_scores_dir is in lowercase
        ld_scores_dir = ld_scores_dir.lower()
        ld_scores_dir = fallExampleDir+"/"+ld_scores_dir
        # Ensure ld_scores_dir has a trailing slash
        if not ld_scores_dir.endswith('/'):
            ld_scores_dir += '/'
        # First command
        command1 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path} --merge-alleles {w_hm3_snplist}  --out {base_name}"
        command12 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path2} --merge-alleles {w_hm3_snplist}  --out {base_name2}"
        
        #command1 = f"python ../munge_sumstats.py --sumstats {sumstats_path} --merge-alleles ../testData/w_hm3.snplist  --out {base_name}"
      
        result1 = subprocess.run( ['bash', '-c', command1], check=True, capture_output=True, text=True)
        result12 = subprocess.run( ['bash', '-c', command12], check=True, capture_output=True, text=True)
      
        print("First command output:", result1.stdout)
        print("First command output:", result12.stdout)
        #print("First command error (if any):", result1.stderr)

        # Second command
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print("Second command ################:",ld_scores_dir)
        command2 = f"cd {fileDir} && python3 {ldsc_script_path} --rg {out_file},{out_file2} --ref-ld-chr {ld_scores_dir} --w-ld-chr {ld_scores_dir} --out {base_name}"
       
        #result2 = subprocess.run( ['bash', '-c', command2], check=True, capture_output=True, text=True)

        try:
            result2 = subprocess.run(['bash', '-c', command2], check=True, capture_output=True, text=True)
            print("Second command output:", result2.stdout)
            separator = "\n--------\n"
            return result1.stdout + separator + result12.stdout + separator + result2.stdout
            #return result2.stdout
        except subprocess.CalledProcessError as e:
            error_details = (e.stderr or e.stdout or e.output or str(e)).strip()
            print(f"An error occurred while running the second command: {command2}")
            print(f"Second command exit code: {e.returncode}")
            print(f"Command stdout: {e.stdout}")
            print(f"Command stderr: {e.stderr}")
            return f"An error occurred while running the second command (exit {e.returncode}): {error_details}"

    except subprocess.CalledProcessError as e:
        error_details = (e.stderr or e.stdout or e.output or str(e)).strip()
        print(f"An error occurred while running the command: {e}")
        print(f"Command stdout: {e.stdout}")
        print(f"Command stderr: {e.stderr}")
        return f"An error occurred while running the command (exit {e.returncode}): {error_details}"


# Example usage
if __name__ == "__main__":
    user_input_sumstats = '../testData/sample/BBJ_HDLC.txt'  # Replace with actual user input
    user_input_ld_scores = '../testData/seu/'  # Replace with actual user input
    #run_herit_command(user_input_sumstats, user_input_ld_scores,False)
    #run_herit_command(user_input_sumstats, user_input_ld_scores,False)