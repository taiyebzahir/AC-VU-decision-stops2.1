import numpy,random

secretion=numpy.array([0,0,0,0,0,1,2])
reception=numpy.array([0,0,0,1,2,0,0])

base=numpy.array( [0,0,1,0,0,0,0])
bound=numpy.array( [1,1,1,1,1,1,1])
SECR=1
ITER=1000
LEAK=0
POP_SIZE=2
cell1, cell2, err, not_complete = 0, 0, 0, 0

m=10000 # number of simulations

trans_mat = numpy.array([[0,-10,0,0,0,0,10], #notch
                         [0,0,0,0,0,1,0], #Delta
                         [0,1,0,0,0,0,0], #basal
                         [0.005,0,0,0,0,0,0], #delta receptor
                         [-10,0,0,0,0,0,0], #notch receptor
                         [0,0,0,0,0,0,0], #ligand_delta
                         [0,0,0,0,0,0,0], #ligand_notch
                         ])



numpy.seterr(all="ignore")

class Stops:
    def __init__(self,mat,pop_size=2,\
         base=None,bound=None,secretion=None,reception=None,\
         secr_amount=1,leak=0,init_env=0):
        self.mat=mat
        self.pop_size=pop_size
        self.secr_amount=secr_amount
        self.leak=leak
        self.init_env=init_env
    
        if base!=None:
            self.base=base
        else: #no base specified
            self.base=numpy.zeros(self.mat.shape[0])

        if bound!=None:
            self.bound=bound
        else: #no bound specified
            self.bound=numpy.ones(self.mat.shape[0])    

        if secretion!=None:
            self.secretion=secretion
        else: #no secretion specified
            self.secretion=numpy.zeros(self.mat.shape[0])

        if reception!=None:
            self.reception=reception
        else: #no reception specified
            self.reception=numpy.zeros(self.mat.shape[0])

        self.init_pop()

    def init_pop(self):
        self.current_pop_size=self.pop_size
        self.pop=numpy.array([list(self.base)]*self.current_pop_size)
        self.env=[]
        for i in range(self.pop_size):
            self.env.append(dict([(x,self.init_env) for x in set(list(self.secretion))-set([0])]))
            #each cell has its own environment and only neighboring cells can receive ligands from this environment

    def can_receive(self,i,j):
        a=numpy.array([\
                       [0,1], 
                       [1,0]
                       ])
        if a[i,j]==1:
            return True
        else:
            return False
    
    def step(self):# gene net change
        tokens_mat=self.pop.dot(self.mat)
        rnd_mat=numpy.random.random(self.pop.shape[0]) # random number for each cell
        sel_mat=numpy.cumsum(abs(tokens_mat),axis=1) # cumulative influence by cell
        sums=numpy.sum(abs(tokens_mat),axis=1) # total influence by cell 
        norm_mat=(numpy.array(sel_mat,dtype=numpy.float32).T/numpy.array(sums)).T # normalized influence by cell
        #print tokens_mat.shape,tokens_mat[:10],rnd_mat.shape,sel_mat,type(tokens_mat),sums.shape
        rnd_mat.resize((self.pop.shape[0],1)) #as a vertical vector
        bool_mat=(norm_mat-rnd_mat)>0 # boolean matrix with values greater than random
        ind_mat=numpy.resize(numpy.array(range(self.pop.shape[1])*self.pop.shape[0])+1,self.pop.shape)
        # matrix of indices
        sel_arr=numpy.select(list(bool_mat.transpose()),list(ind_mat.transpose()))-1
        # the index of the first value greater than random (-1 if no such value)
        dir_arr=numpy.select(list(bool_mat.transpose()),list(numpy.array(tokens_mat).transpose()))
        #
        for i,(s,d) in enumerate(zip(sel_arr,dir_arr)):
            if s>=0:
                self.pop[i,s]=max(0,min(self.bound[s],self.pop[i,s]+(d/abs(d)))) 
             
        #secretion 
        for i,row in enumerate(self.pop):
            secr=row*self.secretion
            self.pop[i]-=secr/secr 
            for j in secr:
                if j!=0:
                    self.env[i][j]=min(1,(self.env[i][j]+self.secr_amount))
 

    #reception
        for i in range(len(self.env)):
            for lig,num in self.env[i].items():
                for lig_num in range(num):
                    index=list(self.reception).index(lig)
                    for j,row in enumerate(self.pop):
                        if self.can_receive(i,j):
                            num-=1
                            self.pop[j,index]=min(self.pop[j,index]+1,self.bound[index])
            

    def sim(self,steps=10000):
        global cell1, cell2, err, not_complete
        for i in range(steps):
            self.step()
        if (self.pop[0,0]==1 and self.pop[1,0]==0) :
            cell1+=1
        elif (self.pop[1,0]==1 and self.pop[0,0]==0) :
            cell2+=1
        elif (self.pop[1,0]==1 and self.pop[0,0]==1) :
            err+=1
        else:
            not_complete+=1

for i in range(m):
    print i
    STOP=Stops(trans_mat,POP_SIZE,base,bound,secretion,reception,secr_amount=SECR,leak=LEAK)
    STOP.sim(ITER)
print "right cell=", cell1, "left cell =", cell2, "error=", err, "not complete=", not_complete
