import numpy as np
import gurobipy
import pandas as pd
from gurobipy import quicksum as sum
import matplotlib.pyplot as plt
def fun_dzq(value_P_wt2h_2,value_P_pv2h_2,rho_wt,rho_pv,lambda_wt,lambda_pv):
    N=range(0,24)
    P_wt2h =np.array( [0, 0, 0, 0, 0, 0, 0.000118965046880248, 1047.26731848920, 1047.26724108527, 1047.26712789578,
              1047.26702135156, 1047.26691403564, 1047.26680723930, 1047.26670219614, 824.427467809682,
              786.259531081652, 1047.26638935748, 1047.26621375800, 1047.26609919524, 1047.26624366009,
              1047.26610449827, 1047.26598842897, 0, 0])
    P_pv2h =np.array([0, 0, 0, 0, 0, 0, 0, 47.2692475009289, 47.2691745816221, 47.2691054587244, 47.2690389348969,
              47.2689731137088, 47.2689067230860, 47.2688424477525, 47.2701311975193, 47.2703042915214,
              47.2686462019993, 47.2685518537297, 47.2684560245753, 0, 0, 0, 0, 0])

    MODEL=gurobipy.Model()
    P_wt_2=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_wt_2',lb=0.34)
    P_pv_2=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv_2',lb=0.40)
    W1=MODEL.addVars(6,vtype=gurobipy.GRB.CONTINUOUS,name='W1')

    MODEL.update()

    MODEL.setObjective(-W1[5]+W1[1]+W1[2]+W1[3]+W1[4])
    MODEL.addConstr(W1[1] == sum(lambda_wt[i] * (P_wt_2[i]-value_P_wt2h_2[i]) for i in N ))
    MODEL.addConstr(W1[2] == sum(lambda_pv[i] * (P_pv_2[i] - value_P_pv2h_2[i] ) for i in N))
    MODEL.addConstr(W1[3] == sum(rho_wt / 2 * (P_wt_2[i] - value_P_wt2h_2[i]) **2 for i in N ))
    MODEL.addConstr(W1[4] == sum(rho_pv / 2 * (P_pv_2[i]- value_P_pv2h_2[i]) for i in N))
    MODEL.addGenConstrLog(W1[0],W1[5])
    MODEL.addConstr(W1[0] == -4903+15743-sum(P_wt_2[i] * P_wt2h[i]+P_pv_2[i] * P_pv2h[i] for i in N))
    #MODEL.addConstrs(-4903+15743-sum(P_wt_2[i] * P_wt2h[i]+P_pv_2[i] * P_pv2h[i] for i in N) >=0)

    #MODEL.Params.NonConvex = 2
    MODEL.optimize()

    obj1=MODEL.objVal   #优化函数值
    rP_pv_2=[P_pv_2[i].x for i in  N]
    rP_wt_2=[P_wt_2[i].x for i in N]
    dict={'obj':obj1,'rP_pv_2':rP_pv_2,'rP_wt_2':rP_wt_2}
    return dict


def fun_fd(P_wt2h_2,rho_wt,lambda_wt):
    P_wt2h =np.array( [0, 0, 0, 0, 0, 0, 0.000118965046880248, 1047.26731848920, 1047.26724108527, 1047.26712789578,
              1047.26702135156, 1047.26691403564, 1047.26680723930, 1047.26670219614, 824.427467809682,
              786.259531081652, 1047.26638935748, 1047.26621375800, 1047.26609919524, 1047.26624366009,
              1047.26610449827, 1047.26598842897, 0, 0])
    P_pv2h =np.array( [0, 0, 0, 0, 0, 0, 0, 47.2692475009289, 47.2691745816221, 47.2691054587244, 47.2690389348969,
              47.2689731137088, 47.2689067230860, 47.2688424477525, 47.2701311975193, 47.2703042915214,
              47.2686462019993, 47.2685518537297, 47.2684560245753, 0, 0, 0, 0, 0])

    MOEEL=gurobipy.Model()

    P_wt_1=MOEEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_wt_1',lb=0.34)
    W=MOEEL.addVars(3,vtype=gurobipy.GRB.CONTINUOUS,name='W')
    MOEEL.update()

    MOEEL.setObjective(-np.log(9757+W[0]-15553)+W[1]+W[2])
    MOEEL.addConstr(W[0] == sum(P_wt2h[i] * P_pv2h[i] for i in range(24)))
    MOEEL.addConstr(W[1] == lambda_wt[i] * (P_wt2h_2[i] - P_wt_1[i])  for i in range(24))
    MOEEL.addConstr(W[2] == rho_wt * (P_wt2h_2[i]- P_wt_1[i]) ** 2 for i in range(24))
    MOEEL.addConstrs(W[0]+9757-15553 >=0)


    MOEEL.Params.NonConvex = 2
    MOEEL.optimize()
    obj=MOEEL.objVal
    ValuePwt2h=[P_wt2h[i].x for i in range(24)]
    dict={'obj':obj,'VP':ValuePwt2h}
    return dict


def fun_gf(P_pv2h_2,rho_pv,lambda_pv):
    P_wt2h =np.array( [0, 0, 0, 0, 0, 0, 0.000118965046880248, 1047.26731848920, 1047.26724108527, 1047.26712789578,
              1047.26702135156, 1047.26691403564, 1047.26680723930, 1047.26670219614, 824.427467809682,
              786.259531081652, 1047.26638935748, 1047.26621375800, 1047.26609919524, 1047.26624366009,
              1047.26610449827, 1047.26598842897, 0, 0])
    P_pv2h =np.array( [0, 0, 0, 0, 0, 0, 0, 47.2692475009289, 47.2691745816221, 47.2691054587244, 47.2690389348969,
              47.2689731137088, 47.2689067230860, 47.2688424477525, 47.2701311975193, 47.2703042915214,
              47.2686462019993, 47.2685518537297, 47.2684560245753, 0, 0, 0, 0, 0])
    MODEL=gurobipy.Model()
    P_pv_1=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv_1',lb=0.4)
    W = MODEL.addVars(5,vtype=gurobipy.GRB.CONTINUOUS,name='W')
    MODEL.update()

    MODEL.setObjective(-W[4]+W[1]+W[2])

    MODEL.addConstr(W[0] == sum(P_pv_1[i] * P_pv2h[i] for i in range(24)))
    MODEL.addConstr(W[1] == sum(lambda_pv[i] * (P_pv2h[i]-P_pv_1[i]) for i in range(24)))
    MODEL.addConstr(W[2] ==sum( rho_pv * (P_pv2h[i]-P_pv_1[i]) ** 2 for i in range(24)))
    MODEL.addConstr(W[3] == 5698 + W[0] - 5932)
    MODEL.addGenConstrLog(W[3], W[4])

    MODEL.Params.NonConvex = 2
    MODEL.optimize()
    obj=MODEL.objVal
    pvh=[P_pv_1[i].x for i in range(24)]
    dict={'obj':obj,'pvh':pvh}
    return dict


def ca():
    P_wt2h =np.array( [0, 0, 0, 0, 0, 0, 0.000118965046880248, 1047.26731848920, 1047.26724108527, 1047.26712789578,
              1047.26702135156, 1047.26691403564, 1047.26680723930, 1047.26670219614, 824.427467809682,
              786.259531081652, 1047.26638935748, 1047.26621375800, 1047.26609919524, 1047.26624366009,
              1047.26610449827, 1047.26598842897, 0, 0])
    P_pv2h =np.array( [0, 0, 0, 0, 0, 0, 0, 47.2692475009289, 47.2691745816221, 47.2691054587244, 47.2690389348969,
              47.2689731137088, 47.2689067230860, 47.2688424477525, 47.2701311975193, 47.2703042915214,
              47.2686462019993, 47.2685518537297, 47.2684560245753, 0, 0, 0, 0, 0])

    ## ADMM迭代参数设置
    rho_wt=1 #惩罚因子
    rho_pv=1 #惩罚因子
    lambda_wt=0*np.ones(24) #风主体拉格朗日乘子
    lambda_pv=0*np.ones(24) #光主体拉格朗日乘子
    maxIter=50 #最大迭代次数
    tolerant=1e-4 #收敛精度
    iter=1 #迭代次数
    #Ben_Store=[] #历史目标函数
    toler1=[] #残差1，风电主体
    toler2=[] #残差2，光伏主体
    P_pv2h_2=np.zeros([maxIter+1,24]) #电制氢主体向光伏主体的期望购电量
    P_wt2h_2=np.zeros([maxIter+1,24]) #电制氢主体向风电主体的期望购电量
    value_P_wt2h_2=np.zeros([maxIter+1,24]) #风电主体向电制氢主体的期望售电量
    value_P_pv2h_2=np.zeros([maxIter+1,24]) #光伏主体向电制氢主体的期望售电量
    Ben_Store1=[]
    Ben_Store2=[]
    Ben_Store3=[]
    P_wt_2=[]
    P_pv_2=[]
    #dict1 return obj Ppv Pwt   dict2 return obj pvh  dict3 return obj VP
    for i in range(maxIter):
        if i == 0:
            dict1=fun_dzq(value_P_wt2h_2[0,:],value_P_pv2h_2[0,:],rho_wt,rho_pv,lambda_wt,lambda_pv)
            Ben_Store1.append(dict1['obj'])
            P_wt_2.append(dict1['rP_wt_2'])
            P_pv_2.append(dict1['rP_pv_2'])

            dict2=fun_gf(dict1['rP_pv_2'],rho_pv,lambda_pv)
            Ben_Store2.append(dict2['obj'])

            dict3=fun_fd(dict1['rP_wt_2'],rho_wt,lambda_wt)
            Ben_Store3.append(dict3['obj'])

            lambda_wt=lambda_wt+rho_wt * (np.array(dict1['Pwt'])-np.array(dict3['VP']))
            lambda_pv=lambda_pv+rho_pv * (np.array(dict1['Ppv'])-np.array(dict2['pvh']))
            toler1.append(np.dot(dict1['rP_wt_2'],dict3['VP']) **2 )
            toler2.append(np.dot(dict1['rP_pv_2'],dict2['pvh']) **2 )
        else:
            dict1=fun_dzq(value_P_wt2h_2[i,:],value_P_pv2h_2[i,:],rho_wt,rho_pv,lambda_wt,lambda_pv)
            Ben_Store1.append(dict1['obj'])
            P_wt_2.append(dict1['rP_wt_2'])
            P_pv_2.append(dict1['rP_pv_2'])

            dict2=fun_gf(dict1['rP_pv_2'],rho_pv,lambda_pv)
            Ben_Store2.append(dict2['obj'])

            dict3=fun_fd(dict1['rP_wt_2'],rho_wt,lambda_wt)
            Ben_Store3.append(dict3['obj'])
            lambda_wt = lambda_wt + rho_wt * (np.array(dict1['r_P_wt_2']) - np.array(dict3['VP']))
            lambda_pv = lambda_pv + rho_pv * (np.array(dict1['rP_pv_2']) - np.array(dict2['pvh']))
            toler1.append(np.dot(dict1['rP_wt_2'],dict3['VP']) **2 )
            toler2.append(np.dot(dict1['rP_pv_2'],dict2['pvh']) **2 )
            if (toler1[i] > tolerant) & (toler2[i] > tolerant):
                break





    plt.figure(1)
    plt.plot(toler1)
    plt.figure(2)
    plt.plot(toler2)
    plt.figure(3)
    plt.plot(P_wt_2)
    plt.figure(4)
    plt.plot(P_pv_2)
    plt.show()

    return [[Ben_Store1],[Ben_Store3],[Ben_Store2]]


if __name__=='__main__':
    data=ca()