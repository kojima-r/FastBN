import pandas as pd
import os

os.makedirs("data_organ",exist_ok=True)
os.makedirs("data_all",exist_ok=True)

df_=pd.read_csv("all_disc.txt",sep="\t")
df=df_.drop(columns=["@name", "tissue"])

df.iloc[:,:10].to_csv("data_all/all_disc10.tsv",index=False,sep="\t")
df.iloc[:,:].to_csv("data_all/all_disc.tsv",index=False,sep="\t")

for t in df_["tissue"].unique():
    df=df_[df_["tissue"]==t]
    df=df.drop(columns=["@name", "tissue"])
    df.iloc[:,:].to_csv("data_organ/"+t+".tsv",index=False,sep="\t")
