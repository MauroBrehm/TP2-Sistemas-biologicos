# Sistema de ED para el modelo de VIH
# dCD4/dt​=− aCD4− bCD4⋅V+ aCD4N​
# dCD8/dt​=− cCD8+ dCD4 . V+ cCD8N​
# dV/dt​= eCD4⋅V− fCD8⋅V− U(t)​
#ademas,
# U(t) = gV(t), g = e . CD4N 
def derivadas (CD4,CD8,V,parametros, droga=False):
    
    a= parametros["a"]
    b= parametros["b"]
    c= parametros["c"]
    d= parametros["d"]
    e= parametros["e"]
    f= parametros["f"]
    CD4N= parametros["CD4N"]
    CD8N= parametros["CD8N"]
    g= parametros["g"]

    dcD4= -a*CD4 - b*CD4*V + a*CD4N 
    dcD8= -c*CD8 + d*CD8*V + c*CD8N
        
    if droga:
        dV= e*CD4*V - f*CD8*V - g*V  #simula inyeccion de la droga por funcion U(t)
    else:
        dV= e*CD4*V - f*CD8*V


    return [dcD4, dcD8, dV]