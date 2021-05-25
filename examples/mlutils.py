# -*- coding: utf-8 -*-
"""

@author: Li-Cheng Xu
"""
import numpy as np
from rdkit import Chem
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error,r2_score
from sklearn.model_selection import train_test_split

def shuffle_index(array,random_state=None):
    np.random.seed(random_state)
    index = list(range(len(array)))
    np.random.shuffle(index)
    return index

def select_exp_set(re_smi,metals,tag,target_smi,target_metal=None,rt=False,temp=None,size=10,random_state=None):
    test_index = []
    if random_state != None:
        np.random.seed(random_state)
    tag_distrib_dict = {0.1:[],0.2:[],0.3:[],0.4:[],
                       0.5:[],0.6:[],0.7:[],0.8:[],
                       0.9:[],1.0:[]}
    shuffle_idx = list(range(len(re_smi)))
    np.random.shuffle(shuffle_idx)
    
    for i in shuffle_idx:
        if re_smi[i] == target_smi:
            if target_metal != None and metals[i] == target_metal:
                if rt and temp[i] >= 20 and temp[i] <= 30:
                    if tag[i] <= 0.1 and len(tag_distrib_dict[0.1]) < size:
                        tag_distrib_dict[0.1].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.2 and tag[i] > 0.1 and len(tag_distrib_dict[0.2]) < size:
                        tag_distrib_dict[0.2].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.3 and tag[i] > 0.2 and len(tag_distrib_dict[0.3]) < size:
                        tag_distrib_dict[0.3].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.4 and tag[i] > 0.3 and len(tag_distrib_dict[0.4]) < size:
                        tag_distrib_dict[0.4].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.5 and tag[i] > 0.4 and len(tag_distrib_dict[0.5]) < size:
                        tag_distrib_dict[0.5].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.6 and tag[i] > 0.5 and len(tag_distrib_dict[0.6]) < size:
                        tag_distrib_dict[0.6].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.7 and tag[i] > 0.6 and len(tag_distrib_dict[0.7]) < size:
                        tag_distrib_dict[0.7].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.8 and tag[i] > 0.7 and len(tag_distrib_dict[0.8]) < size:
                        tag_distrib_dict[0.8].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.9 and tag[i] > 0.8 and len(tag_distrib_dict[0.9]) < size:
                        tag_distrib_dict[0.9].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 1.0 and tag[i] > 0.9 and len(tag_distrib_dict[1.0]) < size:
                        tag_distrib_dict[1.0].append(tag[i])
                        test_index.append(i)
                elif rt == False:
                    if tag[i] <= 0.1 and len(tag_distrib_dict[0.1]) < size:
                        tag_distrib_dict[0.1].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.2 and tag[i] > 0.1 and len(tag_distrib_dict[0.2]) < size:
                        tag_distrib_dict[0.2].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.3 and tag[i] > 0.2 and len(tag_distrib_dict[0.3]) < size:
                        tag_distrib_dict[0.3].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.4 and tag[i] > 0.3 and len(tag_distrib_dict[0.4]) < size:
                        tag_distrib_dict[0.4].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.5 and tag[i] > 0.4 and len(tag_distrib_dict[0.5]) < size:
                        tag_distrib_dict[0.5].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.6 and tag[i] > 0.5 and len(tag_distrib_dict[0.6]) < size:
                        tag_distrib_dict[0.6].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.7 and tag[i] > 0.6 and len(tag_distrib_dict[0.7]) < size:
                        tag_distrib_dict[0.7].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.8 and tag[i] > 0.7 and len(tag_distrib_dict[0.8]) < size:
                        tag_distrib_dict[0.8].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.9 and tag[i] > 0.8 and len(tag_distrib_dict[0.9]) < size:
                        tag_distrib_dict[0.9].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 1.0 and tag[i] > 0.9 and len(tag_distrib_dict[1.0]) < size:
                        tag_distrib_dict[1.0].append(tag[i])
                        test_index.append(i)
            elif target_metal == None:
                if rt and temp[i] >= 20 and temp[i] <= 30:
                    if tag[i] <= 0.1 and len(tag_distrib_dict[0.1]) < size:
                        tag_distrib_dict[0.1].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.2 and tag[i] > 0.1 and len(tag_distrib_dict[0.2]) < size:
                        tag_distrib_dict[0.2].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.3 and tag[i] > 0.2 and len(tag_distrib_dict[0.3]) < size:
                        tag_distrib_dict[0.3].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.4 and tag[i] > 0.3 and len(tag_distrib_dict[0.4]) < size:
                        tag_distrib_dict[0.4].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.5 and tag[i] > 0.4 and len(tag_distrib_dict[0.5]) < size:
                        tag_distrib_dict[0.5].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.6 and tag[i] > 0.5 and len(tag_distrib_dict[0.6]) < size:
                        tag_distrib_dict[0.6].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.7 and tag[i] > 0.6 and len(tag_distrib_dict[0.7]) < size:
                        tag_distrib_dict[0.7].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.8 and tag[i] > 0.7 and len(tag_distrib_dict[0.8]) < size:
                        tag_distrib_dict[0.8].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.9 and tag[i] > 0.8 and len(tag_distrib_dict[0.9]) < size:
                        tag_distrib_dict[0.9].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 1.0 and tag[i] > 0.9 and len(tag_distrib_dict[1.0]) < size:
                        tag_distrib_dict[1.0].append(tag[i])
                        test_index.append(i)
                elif rt == False:
                    if tag[i] <= 0.1 and len(tag_distrib_dict[0.1]) < size:
                        tag_distrib_dict[0.1].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.2 and tag[i] > 0.1 and len(tag_distrib_dict[0.2]) < size:
                        tag_distrib_dict[0.2].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.3 and tag[i] > 0.2 and len(tag_distrib_dict[0.3]) < size:
                        tag_distrib_dict[0.3].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.4 and tag[i] > 0.3 and len(tag_distrib_dict[0.4]) < size:
                        tag_distrib_dict[0.4].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.5 and tag[i] > 0.4 and len(tag_distrib_dict[0.5]) < size:
                        tag_distrib_dict[0.5].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.6 and tag[i] > 0.5 and len(tag_distrib_dict[0.6]) < size:
                        tag_distrib_dict[0.6].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.7 and tag[i] > 0.6 and len(tag_distrib_dict[0.7]) < size:
                        tag_distrib_dict[0.7].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.8 and tag[i] > 0.7 and len(tag_distrib_dict[0.8]) < size:
                        tag_distrib_dict[0.8].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 0.9 and tag[i] > 0.8 and len(tag_distrib_dict[0.9]) < size:
                        tag_distrib_dict[0.9].append(tag[i])
                        test_index.append(i)
                    elif tag[i] <= 1.0 and tag[i] > 0.9 and len(tag_distrib_dict[1.0]) < size:
                        tag_distrib_dict[1.0].append(tag[i])
                        test_index.append(i)
    print('target experiment set size: %d'%len(test_index))
    return test_index

def select_related_set(re_smi,metals,exclude_smi,target_metal,related_smi_set_1,related_smi_set_2):
    related_train_index_1 = []
    related_train_index_2 = []
    for i in list(range(len(re_smi))):
        tmp_smi = re_smi[i]
        if tmp_smi == exclude_smi or metals[i] != target_metal:
            continue
        flag_2 = 0
        
        for tmp_related_2 in related_smi_set_2:
            if Chem.MolFromSmiles(tmp_smi).HasSubstructMatch(Chem.MolFromSmiles(tmp_related_2)):
                flag_2 = 1
        if flag_2 == 0:
            for tmp_related_1 in related_smi_set_1:
                if Chem.MolFromSmiles(tmp_smi).HasSubstructMatch(Chem.MolFromSmiles(tmp_related_1)):
                    related_train_index_1.append(i)
                    break
        elif flag_2 == 1:
            related_train_index_2.append(i)
            
    
    related_train_index_1 = list(set(related_train_index_1))
    related_train_index_2 = list(set(related_train_index_2))
    print('related set 1 size: %d, related set 2 size: %d'%(len(related_train_index_1),len(related_train_index_2)))
    return related_train_index_1,related_train_index_2

class small_sample_learning():
    def __init__(self,related_index_1,related_index_2,test_index,test_size=0.5,split_seed=None):
        np.random.seed(split_seed)
        test_shuffle_index = list(range(len(test_index)))
        np.random.shuffle(test_shuffle_index)
        test_index_shuffle = np.array(test_index)[test_shuffle_index]
        self.related_index_1 = related_index_1
        self.related_index_2 = related_index_2
        self.test_index_shuffle_1 = test_index_shuffle[:int(len(test_index_shuffle)*test_size)]
        self.test_index_shuffle_2 = test_index_shuffle[int(len(test_index_shuffle)*test_size):]
    def delta_learning(self,react_desc,tag,model_ensemble=[],tag_scale=1,n_jobs=1):
        assert len(model_ensemble) == 3, 'model_ensemble should contain 3 models'
        model_1 = model_ensemble[0]
        delta_model_2 = model_ensemble[1]
        delta_model_3 = model_ensemble[2]
        related_train_x_1,related_train_y_1 = react_desc[self.related_index_1],tag[self.related_index_1]
        related_train_x_2,related_train_y_2 = react_desc[self.related_index_2],tag[self.related_index_2]
        append_x,append_y = react_desc[self.test_index_shuffle_1],tag[self.test_index_shuffle_1]
        external_x,external_y = react_desc[self.test_index_shuffle_2],tag[self.test_index_shuffle_2]
        print('delta model is training...')
        model_1.fit(related_train_x_1,related_train_y_1)
        external_y_pred_1 = model_1.predict(external_x)

        pred_related_2 = model_1.predict(related_train_x_2)
        delta_related_2 = related_train_y_2 - pred_related_2
        delta_model_2.fit(related_train_x_2,delta_related_2)
        external_y_pred_2 = delta_model_2.predict(external_x) + model_1.predict(external_x)

        pred_append_y = model_1.predict(append_x)+delta_model_2.predict(append_x)
        delta_append_y = append_y - pred_append_y
        delta_model_3.fit(append_x,delta_append_y)
        external_y_pred_3 = delta_model_3.predict(external_x) + delta_model_2.predict(external_x) + model_1.predict(external_x)

        mae = mean_absolute_error(external_y,external_y_pred_3)*tag_scale
        r2 = r2_score(external_y,external_y_pred_3)
        print('+++delta learning+++MAE: %.3f, r2_score: %.3f'%(mae,r2))
        return external_y_pred_1,external_y_pred_2,external_y_pred_3,external_y
    def only_exp_learning(self,react_desc,tag,model,tag_scale=1,n_jobs=1):
        append_x,append_y = react_desc[self.test_index_shuffle_1],tag[self.test_index_shuffle_1]
        external_x,external_y = react_desc[self.test_index_shuffle_2],tag[self.test_index_shuffle_2]
        print('model training...')
        model.fit(append_x,append_y)
        external_y_pred = model.predict(external_x)
        mae = mean_absolute_error(external_y,external_y_pred)*tag_scale
        r2 = r2_score(external_y,external_y_pred)
        print('+++training with only experiment set+++MAE: %.3f, r2_score: %.3f'%(mae,r2))
        return external_y_pred,external_y
    def with_related_set_raw(self,react_desc,tag,model,tag_scale=1,n_jobs=1):
        related_train_x_1,related_train_y_1 = react_desc[self.related_index_1],tag[self.related_index_1]
        related_train_x_2,related_train_y_2 = react_desc[self.related_index_2],tag[self.related_index_2]
        append_x,append_y = react_desc[self.test_index_shuffle_1],tag[self.test_index_shuffle_1]
        external_x,external_y = react_desc[self.test_index_shuffle_2],tag[self.test_index_shuffle_2]
        
        train_x = np.concatenate([related_train_x_1,related_train_x_2],axis=0)
        train_y = np.concatenate([related_train_y_1,related_train_y_2],axis=0)
        print('model training...')
        model.fit(train_x,train_y)
        external_y_pred = model.predict(external_x)
        mae = mean_absolute_error(external_y,external_y_pred)*tag_scale
        r2 = r2_score(external_y,external_y_pred)
        print('+++training with experiment and related set+++MAE: %.3f, r2_score: %.3f'%(mae,r2))
        return external_y_pred,external_y
class ML():
    def __init__(self,related_set_1,related_set_2,target_sub_set):
        self.r_1_x,self.r_1_y = related_set_1[0],related_set_1[1]
        self.r_2_x,self.r_2_y = related_set_2[0],related_set_2[1]
        self.t_x,self.t_y = target_sub_set[0],target_sub_set[1]
    def hierarc_learn(self,model_ensemble=[],r_t=10):
        assert len(model_ensemble) == 3, 'model_ensemble should contain 3 models, otherwise please manually modify the "hierarchical_learning.hierarc_learn" module'
        r_1_x,r_1_y = self.r_1_x,self.r_1_y
        r_2_x,r_2_y = self.r_2_x,self.r_2_y
        t_x,t_y = self.t_x,self.t_y
        
        base_model = model_ensemble[0]
        delta_model_1 = model_ensemble[1]
        delta_model_2 = model_ensemble[2]
        
        model_ensemble_list = []
        print('model is training...')
        for r in range(r_t):
            print('++++ No. %2d++++'%r)
            #base_model.random_state = r
            #delta_model_1.random_state = r
            #delta_model_2.random_state = r
            
            base_model.fit(r_1_x,r_1_y)
    
            r_2_p = base_model.predict(r_2_x)
            delta_related_2 = r_2_y - r_2_p
            delta_model_1.fit(r_2_x,delta_related_2)
            t_p = base_model.predict(t_x)+delta_model_1.predict(t_x)
            delta_t_y = t_y - t_p
            delta_model_2.fit(t_x,delta_t_y)
            model_ensemble_list.append([base_model,delta_model_1,delta_model_2])

        return model_ensemble_list
    def naive_multi_set_learn(self,model):
        tot_x,tot_y = np.concatenate([self.r_1_x,self.r_2_x,self.t_x],axis=0),np.concatenate([self.r_1_y,self.r_2_y,self.t_y],axis=0)
        print('model is training...')
        model.fit(tot_x,tot_y)
        return model
    def naive_learn(self,model):
        print('model is training...')
        model.fit(self.t_x,self.t_y)
        return model
class eval_models():
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def eval_hierarchic_models(self,model_ensemble_list,scale=1):
        
        total_r2 = []
        total_mae = []
        total_pred_y = []
        for idx,model_ensemble in enumerate(model_ensemble_list):
        
            pred_y = model_ensemble[0].predict(self.x)+\
                            model_ensemble[1].predict(self.x)+model_ensemble[2].predict(self.x)
            tmp_r2 = r2_score(self.y,pred_y)
            tmp_mae = mean_absolute_error(self.y,pred_y)*scale
            total_pred_y.append(pred_y)
            total_r2.append(tmp_r2)
            total_mae.append(tmp_mae)
        highest_r2_idx = np.argmax(total_r2)
        
        print('+++hierarchical learning+++MAE: %.3f, r2_score: %.3f, %d'%(total_mae[highest_r2_idx],total_r2[highest_r2_idx],highest_r2_idx))
        return total_pred_y[highest_r2_idx],model_ensemble_list[highest_r2_idx]
    def eval_naive_model(self,model,scale=1):
        pred_y = model.predict(self.x)
        tmp_r2 = r2_score(self.y,pred_y)
        tmp_mae = mean_absolute_error(self.y,pred_y)*scale
        print('MAE: %.3f, r2_score: %.3f'%(tmp_mae,tmp_r2))
        return pred_y
def draw4fig(ext_y_true,ext_y_pred_exp,ext_y_pred_raw,ext_y_pred_3,ext_y_pred_2,ext_y_pred_1,tag_scale=1,figsave_path=None):
    fig = plt.figure(figsize=(10,8))
    label_font_size = 13
    title_fontsize = 15
    ticks_font_size = 12
    plt.subplot(221)
    plt.scatter(ext_y_true*tag_scale,ext_y_pred_exp*tag_scale,c='lightcoral',alpha=0.8)
    plt.xlabel('Observed $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.ylabel('Predict $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.text(0.2,tag_scale,'MAE: %.3f kcal/mol'%mean_absolute_error(ext_y_true*tag_scale,ext_y_pred_exp*tag_scale),fontsize=ticks_font_size)
    plt.text(0.2,tag_scale-0.5,'${R^2}$: %.3f'%r2_score(ext_y_true*tag_scale,ext_y_pred_exp*tag_scale),fontsize=ticks_font_size)
    plt.plot([0,tag_scale+0.2],[0,tag_scale+0.2],color='lightgrey')
    plt.xticks(fontsize=ticks_font_size)
    plt.yticks(fontsize=ticks_font_size)
    plt.title('prediction performance of set A',fontsize=title_fontsize)
    plt.subplot(222)
    plt.scatter(ext_y_true*tag_scale,ext_y_pred_raw*tag_scale,c='lightblue',alpha=0.8)
    plt.xlabel('Observed $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.ylabel('Predict $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.text(0.2,tag_scale,'MAE: %.3f kcal/mol'%mean_absolute_error(ext_y_true*tag_scale,ext_y_pred_raw*tag_scale),fontsize=ticks_font_size)
    plt.text(0.2,tag_scale-0.5,'${R^2}$: %.3f'%r2_score(ext_y_true*tag_scale,ext_y_pred_raw*tag_scale),fontsize=ticks_font_size)
    plt.plot([0,tag_scale+0.2],[0,tag_scale+0.2],color='lightgrey')
    plt.xticks(fontsize=ticks_font_size)
    plt.yticks(fontsize=ticks_font_size)
    plt.title('prediction performance of set B',fontsize=title_fontsize)
    plt.subplot(223)
    plt.scatter(ext_y_true*tag_scale,ext_y_pred_3*tag_scale,c='yellowgreen',alpha=0.8)
    plt.xlabel('Observed $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.ylabel('Predict $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=label_font_size)
    plt.text(0.2,tag_scale,'MAE: %.3f kcal/mol'%mean_absolute_error(ext_y_true*tag_scale,ext_y_pred_3*tag_scale),fontsize=ticks_font_size)
    plt.text(0.2,tag_scale-0.5,'${R^2}$: %.3f'%r2_score(ext_y_true*tag_scale,ext_y_pred_3*tag_scale),fontsize=ticks_font_size)
    plt.plot([0,tag_scale+0.2],[0,tag_scale+0.2],color='lightgrey')
    plt.xticks(fontsize=ticks_font_size)
    plt.yticks(fontsize=ticks_font_size)
    plt.title('prediction performance of set C',fontsize=title_fontsize)

    plt.tight_layout()
    plt.show()
    if figsave_path != None:
        fig.savefig(figsave_path,dpi=400)
        
def train_eval(info_npz,model,test_size=0.1,tag_scale=1,rand_seed=None,example_mode=False):
    desc = info_npz['desc']
    tag = info_npz['tag']
    if example_mode:
        train_idx = info_npz['train_idx']
        test_idx = info_npz['test_idx']
        train_x,test_x,train_y,test_y = desc[train_idx],desc[test_idx],tag[train_idx],tag[test_idx]
    else:
        np.random.seed(rand_seed)
        train_x,test_x,train_y,test_y = train_test_split(desc,tag,test_size=test_size)
    model.fit(train_x,train_y)
    train_pred = model.predict(train_x)
    test_pred = model.predict(test_x)
    train_r2 = r2_score(train_y,train_pred)
    train_mae = mean_absolute_error(train_y,train_pred)
    test_r2 = r2_score(test_y,test_pred)
    test_mae = mean_absolute_error(test_y,test_pred)
    print('train set MAE: %.3f, r2_score: %.3f'%(mean_absolute_error(train_y,train_pred)*tag_scale,r2_score(train_y,train_pred)))
    print('test set MAE: %.3f, r2_score: %.3f'%(mean_absolute_error(test_y,test_pred)*tag_scale,r2_score(test_y,test_pred)))
    return train_y,train_pred,test_y,test_pred

def drawregfig(train_y,train_pred,test_y,test_pred,tag_scale,figsave_path=None):
    fontsize=18
    fig = plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.scatter(train_y*tag_scale,train_pred*tag_scale,c='darkviolet',alpha=0.2)
    plt.plot([0,4.6],[0,4.6],c='deepskyblue')
    plt.xlabel('Observed $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=fontsize)
    plt.xticks(fontsize=fontsize-3)
    plt.ylabel('Predict $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=fontsize)
    plt.yticks(fontsize=fontsize-3)
    plt.text(0.1,4.4,'MAE: %.3f kcal/mol'%(mean_absolute_error(train_y,train_pred)*tag_scale),fontsize=fontsize)
    plt.text(0.1,4.0,'${R^2}$: %.3f'%r2_score(train_y,train_pred),fontsize=fontsize)

    plt.subplot(122)
    plt.scatter(test_y*tag_scale,test_pred*tag_scale,c='forestgreen',alpha=0.5)
    plt.plot([0,4.6],[0,4.6],c='lightcoral')
    plt.xlabel('Observed $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=fontsize)
    plt.xticks(fontsize=fontsize-3)
    plt.ylabel('Predict $\Delta$$\Delta$$\itG$ (kcal/mol)',fontsize=fontsize)
    plt.yticks(fontsize=fontsize-3)
    plt.text(0.1,4.4,'MAE: %.3f kcal/mol'%(mean_absolute_error(test_y,test_pred)*tag_scale),fontsize=fontsize)
    plt.text(0.1,4.0,'${R^2}$: %.3f'%r2_score(test_y,test_pred),fontsize=fontsize)
    plt.tight_layout()
    if figsave_path != None:
        fig.savefig(figsave_path,dpi=400)

