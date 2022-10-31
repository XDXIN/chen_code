import numpy as np
import pandas as pd
import gurobipy

def fun_dzq(value_P_wt2h_2,value_P_pv2h_2,rho_wt,rho_pv,lambda_wt,lambda_pv):
    L_H2 = np.array(
        [13.1632653100000, 12.6530612200000, 13.6734693900000, 10.6122449000000, 16.1224489800000, 27.9591836700000,
         49.7959183700000,
         51.0204081600000, 47.5510204100000, 46.5306122400000, 47.1428571400000, 48.3673469400000, 49.1836734700000,
         51.2244898000000,
         62.4489795900000, 69.5918367300000, 56.3265306100000, 42.6530612200000, 33.0612244900000, 26.1224489800000,
         20.6122449000000,
         17.9591836700000, 16.3265306100000, 14.4897959200000])
    p_hg = np.array(
        [0.3376 for _ in range(7)] + [0.5980 for i in range(4)] + [0.8654 for i in range(3)] + [0.5980 for i in
                                                                                                range(4)] + [0.8654 for
                                                                                                             i in range(
                4)] + [0.3776 for i in range(2)])
    N=range(1,24)

    MODEL=gurobipy.Model()
    """
    #X=MODEL.addVars(Pv,Pw,Pe,Ph,Pc,mc,Phg,Ebat,Pbac,Pbat,vtype=gurobipy.GRB.INTEGER,name='X')
    #Y=MODEL.addVars(Uabc,Urel,vtype=gurobipy.GRB.BINARY,name='Y')

    P_H2 = sdpvar(1, 24) # % 电制氢主体的产氢量
    P_com = sdpvar(1, 24) #; % 氢气压缩机消耗电功率
    m_com = sdpvar(1, 24) #; % 压缩机压缩氢气流量
    P_hg = sdpvar(1, 24) #; % 电制氢主体从大电网购电量
    E_bat = sdpvar(1, 24) # ; % 储电容量
    P_batc = sdpvar(1, 24) #; % 充电功率
    P_batd = sdpvar(1, 24) #; % 放电功率
    U_abs = binvar(1, 24) #; % 储电设备的放电状态位, 取1时为放电, 0为未放电
    U_relea = binvar(1, 24) #; % 储电设备的充电状态位, 取1时为充电, 0为未充电
    """

    P_pv2h_2 =MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_pv2h_2")    #sdpvar(1, 24) #% 光伏主体向氢气主体购电量（电制氢主体所期望的）
    P_wt2h_2 = MODEL.addVars(1, 24, vtype=gurobipy.GRB.CONTINUOUS, name="P_wt2h_2")  # % 风电主体向氢气主体购电量（电制氢主体所期望的）
    P_el=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_el",lb=0,ub=5000)  #% 产氢对应的耗电量
    P_H2 = MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_H2",lb=0,ub=35)
    P_com=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_com",lb=0,ub=3000)
    m_com=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="m_com",lb=0)
    P_hg=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_hg",lb=0)
    E_bat=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="E_bat",lb=200,ub=1800)
    P_batc=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_batc",lb=0,ub=500)
    P_batd=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_batd",lb=0,ub=600)
    U_abs=MODEL.addVars(1,24,vtype=gurobipy.GRB.BINARY,name="U_abs")
    U_res=MODEL.addVars(1,24,vtype=gurobipy.GRB.BINARY,name="U_res")
    #print("变量创建完成")
    #更新变量
    MODEL.update()


    MODEL.setObjective(np.dot(P_hg,p_hg)+sum(np.sum([0.022*P_el,(1.8e-4)*pow(P_batc/1000,2),(1.8e-4)*pow(P_batd/1000,2)],axis=0).tolist())+
                       np.dot(lambda_wt,[P_wt2h_2[i]-value_P_wt2h_2[i] for i in range(24)])+np.dot(lambda_pv,[P_pv2h_2[i]-value_P_pv2h_2[i] for i in range(24)])+
                       rho_pv/2*pow(np.linalg.norm(np.array(P_pv2h_2)-np.array(value_P_pv2h_2)),2)+rho_wt/2*pow(np.linalg.norm(np.array(P_wt2h_2)-np.array(value_P_wt2h_2)),2)
                       )
    #print("目标函数设置完成")
    #氢气约束
    MODEL.addConstr(P_H2[i]-0.019224*P_el[i]>=0 for i in N)
    MODEL.addVars(P_H2[i]-0.019224*P_el[i]<=0 for i in N)
    #功率约束
    MODEL.addConstr(P_com[i] - 0.2932 * m_com[i] >= 0 for i in N)
    MODEL.addConstrs(P_com[i] - 0.2932 * m_com[i] <= 0 for i in N)
    MODEL.addConstrs(P_H2[0]>=25)
    MODEL.addConstrs(P_H2[0]<=25)
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
    MODEL.addConstrs(E_bat[0] >= 500 + 0.95 * P_batc[0] - P_batd[0] / 0.96)
    MODEL.addConstrs(E_bat[0] <= 500 + 0.95 * P_batc[0] - P_batd[0] / 0.96)

    MODEL.addConstrs(E_bat[i] <= E_bat[i - 1] + 0.95 * P_batc[i] - P_batd[i] / 0.96 for i in range(1, 24))
    MODEL.addConstrs(E_bat[i] >= E_bat[i - 1] + 0.95 * P_batc[i] - P_batd[i] / 0.96 for i in range(1, 24))
    MODEL.addConstrs(E_bat[23]<=500)
    MODEL.addConstrs(E_bat[23]>=500)
    MODEL.addConstrs(P_batc[i]<=U_abs[i]*800 for i in N)
    MODEL.addConstrs(P_batd[i]<=U_res[i]*800 for i in N)
    MODEL.addConstrs(U_abs[i]+U_res[i]<=1 for i in N)


    MODEL.optimize()

    obj1=MODEL.objVal   #优化函数值
    Ppv=[P_pv2h_2[i].x for i in  N]
    Pwt=[P_wt2h_2[i].x for i in N]
    dict={'obj':obj1,'Ppv':Ppv,'Pwt':Pwt}

    return dict


def fun_fd(P_wt2h_2,rho_wt,lambda_wt):
    MOEEL=gurobipy.Model()
    P_wt2g=MOEEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt2g")
    P_wt2h=MOEEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt2h",lb=0)
    P_wt=MOEEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name="P_wt")
    Pwind1max = 1000 *np.array( [1.83969465600000, 2.12213740500000, 1.85496183200000, 2.45038167900000, 2.21374045800000,
                        2.13740458000000, 2.09160305300000, 2.50381679400000, 2.01526717600000, 1.83969465600000,
                        2.36641221400000, 1.90076335900000, 2.07633587800000, 1.71755725200000, 0.824427481000000,
                        0.786259542000000, 1.75572519100000, 2.25190839700000, 1.68702290100000, 2.10687022900000,
                        1.94656488500000, 2.01526717600000, 2.29770992400000, 2.04580152700000])

    MOEEL.update()

    MOEEL.setObjective(-sum(0.34*np.array(P_wt2g))+sum(0.008*np.array(P_wt))+sum(3e-5*pow(np.array(P_wt2h),2)+0.01*P_wt2h)+np.dot(lambda_wt,np.array(P_wt2h)-np.array(P_wt2h_2))
                       +rho_wt / 2 *pow(np.linalg.norm(np.array(P_wt2h)-np.array(P_wt2h_2)),2)
    )
    '''
    C = [];
    C = [C,
         0 <= P_wt2h_1 <= P_wt, % 式(6)：售电不等式约束
    0 <= P_wt2g <= P_wt,
    P_wt2h_1 + P_wt2g == P_wt, % 式(7)：售电等式约束
    0 <= P_wt <= Pwind1max, % 补：风电最大出力约束
    ];'''
    MOEEL.addConstrs(P_wt2h[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt2g[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt2h[i]+P_wt2g[i]<=P_wt[i] for i in range(24))
    MOEEL.addConstrs(P_wt[i]<=Pwind1max[i] for i in range(24))

    MOEEL.optimize()
    obj=MOEEL.objValue
    ValuePwt2h=[P_wt2h[i].x for i in range(24)]
    dict={'obj':obj,'VP':ValuePwt2h}
    return dict


def fun_gf(P_pv2h_2,rho_pv,lambda_pv):
    Ppv1max = 1000 * np.arry([0, 0, 0, 0, 0, 0.244274809000000, 0.435114504000000, 0.740458015000000, 1.12977099200000,
                      1.55725190800000, 1.90839694700000, 1.70992366400000, 1.74045801500000, 1.71755725200000,
                      1.56488549600000, 1.07633587800000, 0.740458015000000, 0.419847328000000, 0.167938931000000, 0, 0,
                      0, 0, 0])
    MODEL=gurobipy.Model()
    P_pv2g=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv2g')
    P_pv2h=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv2h')
    P_pv=MODEL.addVars(1,24,vtype=gurobipy.GRB.CONTINUOUS,name='P_pv')

    MODEL.update()
    MODEL.setObjective(-sum(0.4*np.array(P_pv2h))+sum(0.0085*np.array(P_pv))+sum(3e-5*pow(np.array(P_pv2h),2)+0.01*np.array(P_pv2h))
                       +np.dot(lambda_pv,np.array(P_pv2h_2)-np.array(P_pv2h))+rho_pv / 2 *pow(np.linalg.norm(np.array(P_pv2h_2)-np.array(P_pv2h)),2))

    MODEL.addConstrs(P_pv2h[i]<=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2g[i]<=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2h[i]+P_pv2g[i] <=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv2h[i]+P_pv2g[i] >=P_pv[i] for i in range(24))
    MODEL.addConstrs(P_pv[i] <=Ppv1max[i] for i in range(24))

    MODEL.optimize()

    obj=MODEL.objValue
    pvh=[P_pv2h[i].x for i in range(24)]
    dict={'obj':obj,'pvh':pvh}
    return dict



## ADMM迭代参数设置
rho_wt=1e-4 #惩罚因子
rho_pv=1e-4 #惩罚因子
lambda_wt=0*np.ones(24) #风主体拉格朗日乘子
lambda_pv=0*np.ones(24) #光主体拉格朗日乘子
maxIter=50 #最大迭代次数
tolerant=1e-5 #收敛精度
iter=1 #迭代次数
Ben_Store=[] #历史目标函数
toler1=[] #残差1，风电主体
toler2=[] #残差2，光伏主体
P_pv2h_2=np.zeros([maxIter+1,24]) #电制氢主体向光伏主体的期望购电量
P_wt2h_2=np.zeros([maxIter+1,24]) #电制氢主体向风电主体的期望购电量
value_P_wt2h_2=np.zeros([maxIter+1,24]) #风电主体向电制氢主体的期望售电量
value_P_pv2h_2=np.zeros([maxIter+1,24]) #光伏主体向电制氢主体的期望售电量
#dict1 return obj Ppv Pwt   dict2 return obj pvh  dict3 return obj VP

dict1=fun_dzq(value_P_wt2h_2[0,:],value_P_pv2h_2[0,:],rho_wt,rho_pv,lambda_wt,lambda_pv)
dict2=fun_gf(dict1['Ppv'],rho_pv,lambda_pv)
dict3=fun_fd(dict1['Pwt'],rho_wt,lambda_wt)
lambda_wt=lambda_wt+rho_wt * (dict1['Pwt']-dict3['VP'])
lambda_pv=lambda_pv








