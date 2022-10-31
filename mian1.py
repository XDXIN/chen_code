import numpy as np
import pandas as pd
import gurobipy
from gurobipy import quicksum as sum
import matplotlib.pyplot as plt
def fun_dzq(value_P_wt2h_2,value_P_pv2h_2,rho_wt,rho_pv,lambda_wt,lambda_pv):
    N=range(0,24)
    L_Hh2 = np.array(
        [13.1632653100000, 12.6530612200000, 13.6734693900000, 10.6122449000000, 16.1224489800000, 27.9591836700000,
         49.7959183700000,
         51.0204081600000, 47.5510204100000, 46.5306122400000, 47.1428571400000, 48.3673469400000, 49.1836734700000,
         51.2244898000000,
         62.4489795900000, 69.5918367300000, 56.3265306100000, 42.6530612200000, 33.0612244900000, 26.1224489800000,
         20.6122449000000,
         17.9591836700000, 16.3265306100000, 14.4897959200000])
    L_H2={(i):L_Hh2[i] for i in N}
    p_hhg = np.array(
        [0.3376 for _ in range(7)] + [0.5980 for i in range(4)] + [0.8654 for i in range(3)] + [0.5980 for i in
                                                                                                range(4)] + [0.8654 for
                                                                                                             i in range(
                4)] + [0.3776 for i in range(2)])
    p_hg={(i):p_hhg[i] for i in N}
    MODEL=gurobipy.Model()

    P_pv2h_2 =MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_pv2h_2")    #sdpvar(1, 24) #% 光伏主体向氢气主体购电量（电制氢主体所期望的）
    P_wt2h_2 = MODEL.addVars( 24, vtype=gurobipy.GRB.CONTINUOUS, name="P_wt2h_2")  # % 风电主体向氢气主体购电量（电制氢主体所期望的）
    P_el=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_el",lb=0,ub=5000)  #% 产氢对应的耗电量
    P_H2 = MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_H2",lb=0,ub=35)
    P_com=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_com",lb=0,ub=3000)
    m_com=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="m_com",lb=0)
    P_hg=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_hg",lb=0)
    E_bat=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="E_bat",lb=200,ub=1800)
    P_batc=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_batc",lb=0,ub=500)
    P_batd=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_batd",lb=0,ub=600)
    U_abs=MODEL.addVars(24,vtype=gurobipy.GRB.BINARY,name="U_abs")
    U_res=MODEL.addVars(24,vtype=gurobipy.GRB.BINARY,name="U_res")
    #print("变量创建完成")
    #更新变量
    W1=MODEL.addVars(6,vtype=gurobipy.GRB.CONTINUOUS,name='W1')
    #W2 = MODEL.addVars(vtype = gurobipy.GRB.CONTINUOUS, name='W2')
    #W3 = MODEL.addVars(vtype = gurobipy.GRB.CONTINUOUS, name='W3')
    #W4 = MODEL.addVars(vtype = gurobipy.GRB.CONTINUOUS, name='W4')
    #W5 = MODEL.addVars(vtype = gurobipy.GRB.CONTINUOUS, name='W5')
    #W6 = MODEL.addVars(vtype = gurobipy.GRB.CONTINUOUS, name='W6')

    MODEL.update()
    MODEL.setObjective(sum(W1[i] for i in range(6)))
    MODEL.addConstr(W1[0]==sum(p_hg[i] * P_hg[i] for i in N))
    MODEL.addConstr(W1[1]==sum(0.022 * P_el[i] + 1.8e-4 * ((P_batc[i]/1000)**2+(P_batd[i]/1000)**2) for i in N ) )
    MODEL.addConstr(W1[2]==sum(lambda_wt[i] * (P_wt2h_2[i]-value_P_wt2h_2[i]) for i in N))
    MODEL.addConstr(W1[3]==sum(lambda_pv[i] * (P_pv2h_2[i]-value_P_pv2h_2[i]) for i in N))
    MODEL.addConstr(W1[4]==sum(rho_wt /2 * (P_wt2h_2[i]-value_P_wt2h_2[i])**2 for i in N))
    MODEL.addConstr(W1[5] == sum(rho_pv / 2 * (P_pv2h_2 [i] - value_P_pv2h_2 [i]) ** 2 for i in N))


    #print("目标函数设置完成")
    #氢气约束
    MODEL.addConstrs(P_H2[i]>=0.019224*P_el[i] for i in N)
    MODEL.addConstrs(P_H2[i]<=0.019224*P_el[i] for i in N)
    #功率约束
    MODEL.addConstrs(P_com[i] - 0.2932 * m_com[i] >= 0 for i in N)
    MODEL.addConstrs(P_com[i] - 0.2932 * m_com[i] <= 0 for i in N)
    MODEL.addConstr(P_H2[0]>=25)
    MODEL.addConstr(P_H2[0]<=25)
    #氮气功率平衡
    #P_hg + P_wt2h_2 + P_pv2h_2 + P_batd == P_el + P_batc + P_com, % 式(25): 氢气主体功率平衡
    #0 <= P_hg, % 补充定义: 氢气主体不能向电网售电, 避免光电主体的售电收益被转移到氢气主体侧

    MODEL.addConstrs(P_hg[i]+P_wt2h_2[i]+P_pv2h_2[i]+P_batd[i]<=P_el[i]+P_batc[i]+P_com[i] for i in N)
    MODEL.addConstrs(P_hg[i] + P_wt2h_2[i] + P_pv2h_2[i] + P_batd[i] >= P_el[i] + P_batc[i] + P_com[i] for i in N)
    MODEL.addConstrs(P_el[i]-P_el[i-1]>=-1000 for i in range(1,24))
    MODEL.addConstrs(P_el[i]-P_el[i-1]<=1000 for i in  range(1,24))
    MODEL.addConstrs(P_H2[i]<=P_H2[i-1]+858*(m_com[i]-L_H2[i]) for i in range(1,24))
    MODEL.addConstrs(P_H2[i] >= P_H2[i - 1] + 858 * (m_com[i] - L_H2[i]) for i in range(1, 24))

    #MODEL.addConstrs(E_bat[0] >= 500+0.95*P_batc[0]-P_batd[0]/0.96)
    MODEL.addConstr(E_bat[0] == 500 + 0.95 * P_batc[0] - P_batd[0] / 0.96)
    #MODEL.addConstrs(E_bat[0] <= 500 + 0.95 * P_batc[0] - P_batd[0] / 0.96)

    MODEL.addConstrs(E_bat[i] <= E_bat[i - 1] + 0.95 * P_batc[i] - P_batd[i] / 0.96 for i in range(1, 24))
    MODEL.addConstrs(E_bat[i] >= E_bat[i - 1] + 0.95 * P_batc[i] - P_batd[i] / 0.96 for i in range(1, 24))
    MODEL.addConstr(E_bat[23]==500)
    #MODEL.addConstrs(E_bat[23]>=500)
    MODEL.addConstrs(P_batc[i]<=U_abs[i]*800 for i in N)
    MODEL.addConstrs(P_batd[i]<=U_res[i]*800 for i in N)
    MODEL.addConstrs(U_abs[i]+U_res[i]<=1 for i in N)

    MODEL.Params.NonConvex = 2
    MODEL.optimize()

    obj1=MODEL.objVal   #优化函数值
    Ppv=[P_pv2h_2[i].x for i in  N]
    Pwt=[P_wt2h_2[i].x for i in N]
    dict={'obj':obj1,'Ppv':Ppv,'Pwt':Pwt}

    return dict


def fun_fd(P_wt2h_2,rho_wt,lambda_wt):
    MOEEL=gurobipy.Model()
    P_wt2g=MOEEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt2g")
    P_wt2h=MOEEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt2h",lb=0)
    P_wt=MOEEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt")
    W=MOEEL.addVars(5,vtype=gurobipy.GRB.CONTINUOUS,name='W')
    Pwind1max = 1000 *np.array( [1.83969465600000, 2.12213740500000, 1.85496183200000, 2.45038167900000, 2.21374045800000,
                        2.13740458000000, 2.09160305300000, 2.50381679400000, 2.01526717600000, 1.83969465600000,
                        2.36641221400000, 1.90076335900000, 2.07633587800000, 1.71755725200000, 0.824427481000000,
                        0.786259542000000, 1.75572519100000, 2.25190839700000, 1.68702290100000, 2.10687022900000,
                        1.94656488500000, 2.01526717600000, 2.29770992400000, 2.04580152700000])

    MOEEL.update()

    MOEEL.setObjective(sum(W[i] for i in range(5)))

    MOEEL.addConstr(W[0] == -0.34*sum(P_wt2g))
    MOEEL.addConstr(W[1] == 0.08 * sum(P_wt))
    MOEEL.addConstr(W[2] == sum(3e-5 * ((P_wt[i]) ** 2) +0.01 * P_wt2h[i] for i in range(24)))
    MOEEL.addConstr(W[3] == sum(lambda_wt[i] * (P_wt2h[i]-P_wt2h_2[i]) for i in range(24)))
    MOEEL.addConstr(W[4] == sum(rho_wt /2 * (P_wt2h[i]-P_wt2h_2[i])**2 for i in range(24)))
    MOEEL.addConstrs(P_wt2h[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt2g[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt2h[i]+P_wt2g[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt[i]<=Pwind1max[i] for i in range(24))

    MOEEL.Params.NonConvex = 2
    MOEEL.optimize()
    obj=MOEEL.objVal
    ValuePwt2h=[P_wt2h[i].x for i in range(24)]
    dict={'obj':obj,'VP':ValuePwt2h}
    return dict


def fun_gf(P_pv2h_2,rho_pv,lambda_pv):
    Ppv1max = 1000 * np.array([0, 0, 0, 0, 0, 0.244274809000000, 0.435114504000000, 0.740458015000000, 1.12977099200000,
                      1.55725190800000, 1.90839694700000, 1.70992366400000, 1.74045801500000, 1.71755725200000,
                      1.56488549600000, 1.07633587800000, 0.740458015000000, 0.419847328000000, 0.167938931000000, 0, 0,
                      0, 0, 0])
    MODEL=gurobipy.Model()
    P_pv2g=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv2g')
    P_pv2h=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv2h')
    P_pv=MODEL.addVars(24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv')
    W1=MODEL.addVars(4,vtype=gurobipy.GRB.CONTINUOUS,name='W1')
    #W2 = MODEL.addVars(vtype=gurobipy.GRB.CONTINUOUS, name='W2')
    #W3 = MODEL.addVars(vtype=gurobipy.GRB.CONTINUOUS, name='W3')
    #W4 = MODEL.addVars(vtype=gurobipy.GRB.CONTINUOUS, name='W4')

    MODEL.update()
    MODEL.setObjective(sum(W1[i] for i in range(4)))

    MODEL.addConstr(W1[0]==-0.4*sum(P_pv2h[i] for i in range(24)))
    MODEL.addConstr(W1[1]==0.0085 * sum(P_pv[i] for i in range(24)))
    MODEL.addConstr(W1[2]==sum(3e-5 * (P_pv2h[i])**2+0.01 * P_pv2h[i] for i in range(24)))
    MODEL.addConstr(W1[3]==sum(lambda_pv[i] * (P_pv2h_2[i] - P_pv2h[i] ) for i in range(24)))
    MODEL.addConstrs(P_pv2h[i]<=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2g[i]<=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2h[i]+P_pv2g[i] <=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2h[i]+P_pv2g[i] >=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv[i] <=Ppv1max[i] for i in range(24))

    MODEL.Params.NonConvex = 2
    MODEL.optimize()

    obj=MODEL.objVal
    pvh=[P_pv2h[i].x for i in range(24)]
    dict={'obj':obj,'pvh':pvh}
    return dict


def ca():

    ## ADMM迭代参数设置
    rho_wt=1e-4 #惩罚因子
    rho_pv=1e-4 #惩罚因子
    lambda_wt=0*np.ones(24) #风主体拉格朗日乘子
    lambda_pv=0*np.ones(24) #光主体拉格朗日乘子
    maxIter=50 #最大迭代次数
    tolerant=1e-5 #收敛精度
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
    #dict1 return obj Ppv Pwt   dict2 return obj pvh  dict3 return obj VP
    for i in range(maxIter):
        if i==0:
            dict1=fun_dzq(value_P_wt2h_2[0,:],value_P_pv2h_2[0,:],rho_wt,rho_pv,lambda_wt,lambda_pv)
            Ben_Store1.append(dict1['obj'])

            dict2=fun_gf(dict1['Ppv'],rho_pv,lambda_pv)
            Ben_Store2.append(dict2['obj'])

            dict3=fun_fd(dict1['Pwt'],rho_wt,lambda_wt)
            Ben_Store3.append(dict3['obj'])

            lambda_wt=lambda_wt+rho_wt * (np.array(dict1['Pwt'])-np.array(dict3['VP']))
            lambda_pv=lambda_pv+rho_pv * (np.array(dict1['Ppv'])-np.array(dict2['pvh']))
            toler1.append(np.dot(dict1['Pwt'],dict3['VP']) **2 )
            toler2.append(np.dot(dict1['Ppv'],dict2['pvh']) **2 )
        else:
            dict1=fun_dzq(value_P_wt2h_2[i,:],value_P_pv2h_2[i,:],rho_wt,rho_pv,lambda_wt,lambda_pv)
            Ben_Store1.append(dict1['obj'])

            dict2=fun_gf(dict1['Ppv'],rho_pv,lambda_pv)
            Ben_Store2.append(dict2['obj'])

            dict3=fun_fd(dict1['Pwt'],rho_wt,lambda_wt)
            Ben_Store3.append(dict3['obj'])
            lambda_wt = lambda_wt + rho_wt * (np.array(dict1['Pwt']) - np.array(dict3['VP']))
            lambda_pv = lambda_pv + rho_pv * (np.array(dict1['Ppv']) - np.array(dict2['pvh']))
            toler1.append(np.dot(dict1['Pwt'],dict3['VP']) **2 )
            toler2.append(np.dot(dict1['Ppv'],dict2['pvh']) **2 )
            if (toler1[i] > tolerant) & (toler2[i] > tolerant):
                break






    plt.plot(Ben_Store1)
    plt.plot(Ben_Store2)
    plt.plot(Ben_Store3)
    plt.show()

    return [[Ben_Store1],[Ben_Store3],[Ben_Store2]]


if __name__=='__main__':
    data=ca()