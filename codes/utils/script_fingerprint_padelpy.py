import os
import glob
import pandas as pd
from padelpy import padeldescriptor

def get_padel_fp(cid, smile):
    print("The compound CID : " + cid)
    print("The compound SMILES : " + smile)

    dir_fp = "data/middle/fp_padel/" + cid
    tmp_smi = "data/middle/tmp_" + cid + ".smi"
    FP_list = ['AtomPairs2DCount', # 781 +
         'AtomPairs2D', # 781
         'EState', #80 *
         'CDKextended', #1025 *
         'CDK', #1025 *
         'CDKgraphonly', #1025 *
         'KlekotaRothCount', #4861 * +
         'KlekotaRoth', #4861 * +
         'MACCS', #167
         'PubChem', #882
         'SubstructureCount', #308 +
         'Substructure'] #308

    if not os.path.exists(dir_fp):
    # if len(os.listdir(dir_fp)) != 12:
        os.makedirs(dir_fp,exist_ok=True)

        with open(tmp_smi, "w") as f:
            f.write(smile + "\t" + cid)

        xml_files = glob.glob("/home/lishensuo/PROJECT/Toxicity/data/feature/fingerp/fingerprints_xml/*.xml")
        xml_files.sort()
        fp_xml_dict = dict(zip(FP_list, xml_files))
        fp_meta = pd.read_csv("/home/lishensuo/PROJECT/Toxicity/data/feature/merge/All_feat_meta19793.csv")

        for FP in FP_list :
            print(FP)
            # FP = FP_list[0]
            output_file = "data/middle/fp_padel/" + cid + "/" + FP + '.csv'
            try:
                padeldescriptor(mol_dir=tmp_smi, 
                                d_file=output_file, 
                                descriptortypes= fp_xml_dict[FP],
                                detectaromaticity=True,
                                standardizenitro=True,
                                standardizetautomers=True,
                                removesalt=True,
                                log=False,
                                fingerprints=True)
            except:
                print("==> Failed to code, return None")
                fp_meta_sub = fp_meta[fp_meta.Type.isin([FP])] #查询特定指纹的长度
                output_df = pd.DataFrame([cid] + [None]*fp_meta_sub.shape[0]).transpose() 
                output_df.columns = ["Name"] + list(fp_meta_sub.Name)
                output_df.to_csv("data/middle/fp_padel/" + cid + "/" + FP + '.csv', index=False)
        os.remove(tmp_smi)

    fp_feat = pd.DataFrame({"Name":[int(cid)]})
    for FP in FP_list:
        # FP = FP_list[0]
        feat_tmp = pd.read_csv("data/middle/fp_padel/" + cid + "/" + FP + '.csv')
        fp_feat = pd.merge(fp_feat, feat_tmp)
    return fp_feat