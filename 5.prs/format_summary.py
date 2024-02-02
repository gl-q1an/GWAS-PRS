# DESCRIPTION:This script is used to automatically transform the imported GWAS summary into 
#    those columns we want, and try to calculate those columns we want

import os
import numpy as np
import pandas as pd
import json
import argparse
import re
import datetime
from contextlib import redirect_stdout

Header = "=============================================\n"
Header +="||      Automated GWAS Summary Converter   ||\n"
Header +="=============================================\n"
Header += str(datetime.datetime.now())+"\n"

describe_dict = {
    "args_all":Header,
    "describe":"An Automated GWAS Summary Converter",
    "summary":" The GWAS summary you want to standardize, .txt or .gz",
    "out":"     The file name of standardized summary, .txt or .gz",
    "outform":" Which col you want to choose to output,[SNP,CHR,BP,A1,A2,Frq,P,Z,BETA,SE,OR,N,Nca,Nco] default SNP,CHR,A1,A2,Z,N",
    "n_num":"   If there is no column representing N, you can give a exact value for N, which will add a column N to your output",
    "nca_num":" Similar the explanation of n_num",
    "nco_num":" Similar the explanation of N_num",
    "n_cal":"direct or effect<d>. If Nca Nco N exists at the same time, and if N is needed, direct means taking N directly without calculation, and effected means taking d for calculation. The default is effect4,which means d=4,N=d/(1/Nca+1/Nco). If d=0,N=Nca+Nco."
}

##### List of functions for relevant steps #####
# 1. Write the log file
def printlog(text):
    global logfile
    print(text)
    with open(logfile,"a") as f:
        with redirect_stdout(f):
            print(text)
            
# 2. Identify the args and create the log file
def indentify_args(describe_dict):
    parser = argparse.ArgumentParser(description=describe_dict["describe"])
    parser.add_argument("-s","--summary", required=True, help=describe_dict["summary"])
    parser.add_argument("-o", "--out", required=True, help=describe_dict["out"])
    parser.add_argument("-f","--outform", default="SNP,CHR,A1,A2,Z,N", help=describe_dict["outform"])
    parser.add_argument("--n_num", help=describe_dict["n_num"])
    parser.add_argument("--nca_num", type=int, help=describe_dict["nca_num"])
    parser.add_argument("--nco_num", type=int, help=describe_dict["nco_num"])
    parser.add_argument('--n_cal', default="effect4", help=describe_dict["n_cal"])
    args = parser.parse_args()
    global logfile
    logfile=args.out+".log"
    with open(logfile,"w") as f:
        f.write(describe_dict["args_all"])
    printlog(describe_dict["args_all"])
    if not re.match(r'(direct|effect\d+)', args.n_cal):
        printlog("Erorr:n_cal should be `direct` or `effect<d>`, such as effect4! See -h or --help")
        exit()
    return args

# 3. Read the json file
# Automatically identify column names based on information in json files
def identify_json(jsonfile):
    with open(jsonfile, 'r') as json_file:
        load = json.load(json_file)
        ref_dict = {key: value for key, value in load.items() if "DESCRIPTION" not in key}
        res_dict = {key: None for key in ref_dict}
    return ref_dict,res_dict

# 4. Change strings in json to lowercase to improve matching ability
def trans2low(input_str):
    output_str=""
    for char in input_str:
        if char.isalpha() and char.isupper():
            output_str += char.lower()
        else:
            output_str += char
    return output_str

# 5. Read the summary file
def read_summary(inputsum):
    if inputsum.endswith(".gz"):
        df=pd.read_csv(inputsum, compression='gzip',engine='python', sep=None)
    else:
        df=pd.read_csv(inputsum, engine='python',sep=None)
    num_rows, num_columns = df.shape
    printlog(f"--There are {num_rows} SNPs!!--\n")
    if num_columns < 3:
        printlog("Error:Unreasonable number of columns!!")
        exit()
    return df

# 6. Match to the reference
def auto_match(df,ref_dict,res_dict,args):
    for col in df.columns:
        low_col=trans2low(col)
        for k,v in ref_dict.items():
            if low_col in v:
                if res_dict[k] is None:
                    res_dict[k] = col
                else:
                    printlog(f"**Warning:Col {res_dict[col]} and {k} both match {col}. Choose the former.")
    cond1=res_dict["Nca"] is not None and res_dict["Nco"] is not None
    cond2=args.nca_num is not None and args.nco_num is not None
    if re.match(r'effect\d+', args.n_cal) and (cond1 or cond2):
        res_dict["N"] = None
    return res_dict

# 7. Change the colnames and specify a certain value for N
def rename_col(df, res_dict, args):
    match_dict={}
    for k, v in res_dict.items():
        if v is not None:
            match_dict[v] = k
    df = df.rename(columns=match_dict)
    printlog("------Here are cols matched!------")
    for k,v in match_dict.items():
        printlog(f"    {k}--->{v}   ")
    printlog("----------------------------------")
    
    if args.n_num is not None:
        df["N"] = args.n_num
        res_dict["N"] = "N"
        printlog(f"N is specified as {args.n_num}")
    if args.nca_num is not None:
        df["Nca"] = args.nca_num
        res_dict["Nca"] = "Nca"
        printlog(f"Nca is specified as {args.nca_num}")
    if args.nca_num is not None:
        df["Nco"] = args.nco_num
        res_dict["Nco"] = "Nco"
        printlog(f"Nco is specified as {args.nco_num}")
    return df, res_dict

# 8.Identify which col is missing
def missing_cols(res_dict, outform):
    opt=[]
    outformlist=outform.split(",")
    for k,v in res_dict.items():
        if k in outformlist and v is None:
            opt.append(k)
    return opt

# 9.Important! Calculate the missing col!
def cal_misscol(df, opt, res_dict, ref_dict, args):
    
    if "SNP" in opt:
        printlog("-SNP col is not found!")
        if all(col in df.columns for col in ["CHR", "BP", "A1"]):
            printlog("|-Get is from CHR:BP:A1")
        elif all(col in df.columns for col in ["CHR", "BP"]):
            printlog("|-Get is from CHR:BP")
        elif "chrpos" in df.columns:
            for i in ref_dict["SEP_OF_CHR_POS"]:
                if i in str(df.loc[0, 'chrpos']):
                    if i != ":":
                        df['SNP'] = df['chrpos'].str.replace(i, ":")
                        break
            printlog("|-Get is from chrpos col")
        elif "chrposa1a2" in df.columns:
            for i in ref_dict["SEP_OF_CHR_POS"]:
                if i in str(df.loc[0, 'chrposa1a2']):
                    if i != ":":
                        df['SNP'] = df['chrposa1a2'].str.replace(i, ":")
                        break
            printlog("|-Get is from chrposa1a2 col")
        else:
            printlog("Error:SNP col is not found!!")
            exit()
        
    if "CHR" in opt:
        if "chrpos" in df.columns:
            printlog("-CHR col is not found!")
            for i in ref_dict["SEP_OF_CHR_POS"]:
                if i in str(df.loc[0, 'chrpos']):
                    sep=i
                    break
            printlog("|-Get it from chrpos! The sep is "+repr(sep))
            df['CHR'] = df["chrpos"].str.split(sep, expand=True)[0]
            res_dict['CHR']='CHR'
        elif "chrposa1a2" in df.columns:
            printlog("-CHR col is not found!")
            for i in ref_dict["SEP_OF_CHR_POS"]:
                if i in str(df.loc[0, 'chrposa1a2']):
                    sep=i
                    break
            printlog("|-Get it from chrposa1a2! The sep is "+repr(sep))
            df['CHR'] = df["chrposa1a2"].str.split(sep, expand=True)[0]
            res_dict['CHR']='CHR'
        else:
            printlog("Error:CHR col is not found!!")
            printlog("Currently it is not able to automatically obtain information based on SNPs and references, but this may be added later.")
            exit()
            
    if "BP" in opt:
        if "chrpos" in df.columns:
            printlog("-BP col is not found!")
            printlog("|-Get it from chrpos! The sep is "+repr(sep))
            df['BP'] = df["chrpos"].str.split(sep, expand=True)[1]
            res_dict['BP']='BP'
        elif "chrposa1a2" in df.columns:
            printlog("-BP col is not found!")
            printlog("|-Get it from chrposa1a2! The sep is "+repr(sep))
            df['BP'] = df["chrposa1a2"].str.split(sep, expand=True)[1]
            res_dict['BP']='BP'
        else:
            printlog("Error:BP col is not found!!")
            printlog("Currently it is not able to automatically obtain information based on SNPs and references, but this may be added later.")
            exit()
            
    if "A1" in opt:
        if "chrposa1a2" in df.columns:
            printlog("-A1 col is not found!")
            printlog("|-Get it from chrposa1a2! The sep is "+repr(sep)) 
            df['A1'] = df["chrposa1a2"].str.split(sep, expand=True)[2]
            res_dict['A1']='A1'
        else:
            printlog("Error:A1 is not found!!")
            printlog("Currently it is not able to automatically obtain information based on SNPs and references, but this may be added later.")
            exit()
            
    if "A2" in opt:
        if "chrposa1a2" in df.columns:
            printlog("-A2 col is not found!")
            printlog("|-Get it from chrposa1a2! The sep is "+repr(sep)) 
            df['A2'] = df["chrposa1a2"].str.split(sep, expand=True)[3]
            res_dict['A2']='A2'
        else:
            printlog("Error:A2 is not found!!")
            printlog("Currently it is not able to automatically obtain information based on SNPs and references, but this may be added later.")
            exit() 
    
    if "Z" in opt:
        if "BETA" in df.columns and "SE" in df.columns:
            printlog("-Z col is not found!")
            printlog("|-Get it from BETA and SE, Z=BETA/SE!") 
            df['Z'] = (df['BETA']/df['SE']).round(5)
            res_dict['Z']='Z'
        elif "OR" in df.columns and "SE" in df.columns:
            printlog("-Z col is not found!")
            printlog("|-Get it from BETA and SE, Z=log(OR)/SE!") 
            df['Z'] = (df['OR'].apply(lambda x: np.log(x))/df['SE']).round(5)
            res_dict['Z']='Z'
        else:
            printlog("Error:Z is not found!!")
            printlog("Maybe Z can be calculated by P, it may be added later!")
            printlog("We just compute it by BETA/SE but not p")
            exit()        
    
    if "N" in opt:
        if "Nca" in df.columns and "Nco" in df.columns:
            if args.n_cal == "direct":
                printlog("-Get N from Nco and Nca, effect=0, N=Nco+Nca") 
                df['N'] = df['Nca']+df['Nco']
                res_dict['N']='N'
            else:
                d=int(args.n_cal[-1])
                if d==0:
                    printlog(f"-Get N from Nco and Nca, d={d}, N={d}/(1/Nco+1/Nca)")
                    df['N'] = df['Nca']+df['Nco']
                    df['N'] = df['N'].astype(int)
                else:
                    printlog(f"-Get N from Nco and Nca, d={d}, N={d}/(1/Nco+1/Nca)")
                    df['N'] = d/(1/df['Nca']+1/df['Nco'])
                    df['N'] = df['N'].astype(int)
        else:
            printlog("Error:N is not found!!")
            exit() 
    
    if "BETA" in opt:
        if "OR" in df.columns:
            printlog("-BETA col is not found!")
            printlog("|-Get it from OR! BETA=log(OR)") 
            df['BETA'] = (df['OR'].apply(lambda x: np.log(x))).round(5)
            res_dict['BETA']='BETA'
        elif all(col in df.columns for col in ["Z", "Frq", "N"]):
            printlog("-BETA col is not found!")
            printlog("|-Get it from z/d! d=[2f(1-f)(z^2+N)]^0.5")
            df['d_tmp']=(2*df['Frq']*(1-df['Frq'])*(df['Z']**2+df['N']))**0.5
            df['BETA'] = (df['Z']/df['d_tmp']).round(5)
            res_dict['BETA']='BETA'
        else:
            printlog("Error:BETA is not found!!")
            printlog("Maybe BETA can be calculated, it may be added later!")
            printlog("We just compute it by BETA=log(OR) or by Z|Frq|N")
            exit()
    
    if "SE" in opt:
        if 'd_tmp' in df.columns:
            printlog("-SE col is not found!")
            printlog("|-Get it from 1/d! d=[2f(1-f)(z^2+N)]^0.5")
            df['SE'] = (1/df['d_tmp']).round(5)
            res_dict['SE']='SE'
        elif all(col in df.columns for col in ["Z", "Frq", "N"]):
            printlog("-SE col is not found!")
            printlog("|-Get it from 1/d! d=[2f(1-f)(z^2+N)]^0.5")
            df['d_tmp']=(2*df['Frq']*(1-df['Frq'])*(df['Z']**2+df['N']))**0.5
            df['SE'] = (1/df['d_tmp']).round(5)
            res_dict['SE']='SE'
        else:
            printlog("Error:SE is not found!!")
            printlog("Maybe SE can be calculated, it may be added later!")
            printlog("We just compute it by 1/d")
            exit()
            
    if "OR" in opt:
        if "BETA" in df.columns:
            printlog("-OR col is not found!")
            printlog("|-Get it from BETA! OR=exp(BETA)")         
            df['OR'] = (df['BETA'].apply(lambda x: np.exp(x))).round(5)
            res_dict['OR']='OR'
        else:
            printlog("Error:OR is not found!!")
            printlog("Maybe OR can be computed, it may be added later!")
            printlog("We just compute it by OR=exp(BETA)")
            exit()
                
    if "A1" in df.columns:
        if df.at[0,"A1"].islower:
            df['A1'] = df['A1'].str.upper()
    if "A2" in df.columns:
        if df.at[0,"A1"].islower:
            df['A2'] = df['A2'].str.upper()    
        
    return df, res_dict

def output(df,outfile,outform):
    outformlist=outform.split(",")
    try:
        df=df[outformlist]
    except KeyError:
        for i in outformlist:
            if i not in df.columns:
                outformlist.remove(i)
                printlog(f"-{i} is not found!!")
        df=df[outformlist]

    if "CHR" in outformlist and "BP" in outformlist:
        df=df.sort_values(by=["CHR","BP"])
    if outfile.endswith(".gz"):
        df.to_csv(outfile, index=False, sep='\t', compression='gzip')
    else:
        df.to_csv(outfile, index=False, sep='\t')
    printlog(f"{outfile} is writing!!")
# ===== Function list end =====

def main(describe_dict):
    # read the json file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_file=os.path.join(script_dir,"format_match.json")
    ref_dict,res_dict=identify_json(json_file)
    
    # read the args
    args=indentify_args(describe_dict)
    
    # read the summary and match
    df=read_summary(args.summary)
    res_dict=auto_match(df,ref_dict,res_dict,args)
    df,res_dict=rename_col(df, res_dict, args)
    
    # calculate the missing col
    opt=missing_cols(res_dict, args.outform)
    df,res_dict=cal_misscol(df, opt, res_dict, ref_dict, args)
    output(df, args.out, args.outform)
    printlog("Finished!!"+str(datetime.datetime.now()))
    
if __name__=="__main__":
    main(describe_dict)