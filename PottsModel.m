%% Potts model: generalized Ising, 3 states

close all
clear all

% INITIALIZATION of the system. 

%setting J to 1 aand Kb to 1 the transition temperature is around 2 (?)

Kb=1;
%T=[0.1,0.2,0.4,0.6,0.9,1.1,1.5,1.9,2.5,3];
T=[0.2,0.8,1.5];
J=1;


burn_in=3000000;%iterations to reach equilibrium, empirically set
indep_stat=50000;
total_nbr_iter=burn_in+(indep_stat*5);
beta=1./(Kb*T);
side_dim=80;
grid_dim=side_dim*side_dim;

grid_pos=zeros(grid_dim,2);%list of positions of each grid site
spin_list=zeros(grid_dim,1);%list of spin of each grid site
conf_T=zeros(grid_dim,length(T));
magn_av=zeros(1,length(T));
%%

%implement boundary conditions: done by creating a list that 
%contains the (indices of) neighbors of each point.
    
 a=1;

for i=1:side_dim
    for j=1:side_dim
        grid_pos(a,:)=[j,i];
        spin_list(a)=sign(rand-1/2);
        a=a+1;
    end
end


%nearest neighbors: they will be stored in order right, left, up, down.

nn_list=zeros(grid_dim,4);
nn_list(1,:)=[2,side_dim,grid_dim-side_dim+1,side_dim+1];


for a=2:side_dim-1
    nn_list(a,:)=[a+1,a-1,grid_dim-side_dim+a,side_dim+a];
end
nn_list(side_dim,:)=[1,side_dim-1,grid_dim,side_dim*2];

for a=side_dim+1:grid_dim-side_dim-1
    nn_list(a,:)=[a+1,a-1,a-(side_dim),a+(side_dim)];
end

for i=2:side_dim-1
nn_list(side_dim*(i-1)+1,:)=[side_dim*(i-1)+2,side_dim*i,side_dim*(i-2)+1,side_dim*i+1];
nn_list(side_dim*i,:)=[side_dim*(i-1)+1,side_dim*i-1,side_dim*(i-1),side_dim*(i+1)];
end

for a=grid_dim-side_dim+2:grid_dim-1
    nn_list(a,:)=[a+1,a-1,a-side_dim,a-grid_dim+side_dim];
end

nn_list(grid_dim-side_dim+1,:)=[grid_dim-side_dim+2,grid_dim,grid_dim-(side_dim*2)+1,1];
nn_list(grid_dim,:)=[grid_dim-side_dim+1,grid_dim-1,grid_dim-side_dim,side_dim];    




%% ISING model

%the Hamiltonian should be
%H=-J*sum_<ij>(Si*Sj)-sum_i(hi*Si)
%where <ij> are nearest neighbours and hi is the external field that spin i
%feels. 
%we decide there is no external field, so the only forces are given by spin
%interactions. 


%if the piece is ferromagnetic it will have a spontaneous magnetization
%so H<0 for parallel spins, H>0 for antiparallel spins. 



for temp=1:length(T) %final aim: entropy as function of temperature
%% RANDOM INITIALIZATION OF THE SYSTEM AT EVERY TEMPERATURE CHANGE
    


for a=1:grid_dim
        spin_list(a)=sign(rand-1/2);
end

stat=0;

for iter=1:total_nbr_iter   %run continuously, collect statistics only after burn in

a=1+floor(rand*grid_dim);


diff=spin_list(a)-spin_list(nn_list(a,:));
tool=find(~diff);
neighb1=-J*length(tool);


%propose a change random
nn=rand-1/2;
if spin_list(a)==0
    aprop=sign(nn);
elseif spin_list(a)==1
    if nn<=0
    aprop=-1;
    else
        aprop=0;
    end
    
else 
    if nn<=0
       aprop=1;
    else
        aprop=0;
    end 
    
end

diff=aprop-spin_list(nn_list(a,:));
tool=find(~diff);
neighb2=-J*length(tool);




%E_curr=-J*sum(neighb1);
%E_prop=-J*sum(neighb2);
%deltaE=E_prop-E_curr;
deltaE=neighb2-neighb1;

if deltaE<0
    spin_list(a)=aprop;
else 
    prob_acc=exp(-deltaE*beta(temp));
    q=rand;
    if q<prob_acc
        spin_list(a)=aprop;
    end
end


%OTHER WAY TO COMPUTE IT. IT WORKS. SAME
%E_curr=exp(-beta(temp)*(J*sum(neighb1)));
%E_prop=exp(-beta(temp)*(J*sum(neighb2)));
%En_ratio=E_prop/(E_curr+E_prop);
% if En_ratio<1
%     %the proposed state has lower energy, accept it with a prob En_ratio
%     q=rand;
%     if q>En_ratio %accepted, flip spin!
%         spin_list(a)=-spin_list(a);
%     end
% end

  
 

%% COLLECT STATISTICS

if iter==burn_in+stat*indep_stat
   stat=stat+1;

  magnetization(stat)=sum(spin_list)/grid_dim;
  
  end%after burn in


end%end of total loop

magn_av(temp)=sum(magnetization)/stat;




%% SAVE CONFIGURATION 

conf_T(:,temp)=spin_list;
figure 

  for i=1:grid_dim
      if conf_T(i,temp)==1
           plot(grid_pos(i,1),grid_pos(i,2),'r*');hold on
      elseif conf_T(i,temp)==-1
          plot(grid_pos(i,1),grid_pos(i,2),'b*');hold on
      else 
          plot(grid_pos(i,1),grid_pos(i,2),'g*');hold on
       end
   end



end %(temperature cf)


figure
plot(T,magn_av,'bo');


