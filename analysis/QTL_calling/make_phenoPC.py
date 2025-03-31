from sklearn.decomposition import PCA
import pandas as pd
import argparse

#######################
### Author Gracie Gordon 2021
### Calculate PCs for phenotype of interest
#######################

parser = argparse.ArgumentParser(description="Calculate PCs for phenotype of interest.")
parser.add_argument("-p", "--n_pc", type=int, required=True, help='number of PCs')
parser.add_argument("-e", "--pheno", type=str, required=True, help='phenotype matrix')
parser.add_argument("-o", "--outfile", type=str, required=True, help='outfile')

args = parser.parse_args()

def exprPC(file,pc_n,name):
    pca = PCA(n_components=pc_n)
    #e=pd.read_csv(file,sep="\t")
    e=pd.read_csv(file,sep=" ")
    #print(e.head())
    #etmp=e.drop(columns="Unnamed: 0")
    #sample=e["Unnamed: 0"]
    etmp=e
    principalComponents = pca.fit_transform(etmp.T)
    c=list(range(1,(pc_n+1)))
    cols=[str("phenPC"+str(i)) for i in c]
    principalDf = pd.DataFrame(data = principalComponents
             , columns = cols)

    #principalDf.index=e.columns[1:]
    principalDf.index=e.columns
    out=principalDf.T
    #print(out)
    out.to_csv(name,sep="\t",index=True)

    # save variance explained
    print('explained var')
    #print(pca.explained_variance_ratio_.cumsum())
    #print(pca.explained_variance_ratio_)
    var_exp=pd.concat([pd.DataFrame(cols),pd.DataFrame(pca.explained_variance_ratio_)],axis=1)
    var_exp.columns=['PC','variance_explained']
    var_exp.to_csv(str(name+'_vexp.txt'),sep="\t",index=False)
    #var_exp= pd.DataFrame(data=pca.explained_variance_ratio_.cumsum(),columns=cols)
    #var_exp.index=e.columns
    #print(var_exp)

print('npc: '+str(args.n_pc))
print('pheno: '+str(args.pheno))
print('outfile: '+str(args.outfile))
exprPC(args.pheno,args.n_pc,args.outfile)
