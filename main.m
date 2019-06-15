%%Project Information Theory
%%Laura Masaracchia

%% 2-D Spin system 

close all
clear all

% INITIALIZATION of the system. 

%setting J to 1 aand Kb to 1 the transition temperature is around 2 (?)

Kb=1;
T=[0.2,0.8,1.2,1.7,2.1,2.5,2.9,3.5,5,10];
J=1;


burn_in=3000000;%iterations to reach equilibrium, empirically set
indep_stat=50000;
total_nbr_iter=burn_in+(indep_stat*20);
beta=1./(Kb*T);
side_dim=50;
grid_dim=side_dim*side_dim;

grid_pos=zeros(grid_dim,2);%list of positions of each grid site
spin_list=zeros(grid_dim,1);%list of spin of each grid site
conf_T=zeros(grid_dim,length(T));

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
%H=-J*sum_<ij>delta(Si*Sj)-sum_i(hi*Si)
%any even couple gives a + contribution, any uneven couple gives a zero
%so the min can be 5 spins equal, -4J
%
%where <ij> are nearest neighbours and hi is the external field that spin i
%feels. 
%we decide there is no external field, so the only forces are given by spin
%interactions. 


%if the piece is ferromagnetic it will have a spontaneous magnetization
%so H<0 for parallel spins, H>0 for antiparallel spins. 



for temp=1:length(T) %final aim: entropy as function of temperature
%% RANDOM INITIALIZATION OF THE SYSTEM AT EVERY TEMPERATURE CHANGE
   
%three states

for a=1:grid_dim
    q=rand*3;
    if q<=1
        spin_list(a)=-1;
    elseif q<=2
        spin_list(a)=0;
    else
        spin_list(a)=+1;
    end
end

stat=0;

for iter=1:total_nbr_iter   %run continuously, collect statistics only after burn in

a=1+floor(rand*grid_dim);

neighb1=spin_list(nn_list(a,:)).*spin_list(a);
neighb2=spin_list(nn_list(a,:)).*(-spin_list(a));


E_curr=-J*sum(neighb1);
E_prop=-J*sum(neighb2);
deltaE=E_prop-E_curr;
if deltaE<0
    spin_list(a)=-spin_list(a);
else 
    prob_acc=exp(-deltaE*beta(temp));
    q=rand;
    if q<prob_acc
        spin_list(a)=-spin_list(a);
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



    nbr_block1=zeros(2,4);%nbr_block1(1,1) will contain block1 A up
                      %nbr_block1(2,1) will contain block1 A down
                      % . . . 
                      %nbr_block1(1,4) will contain block1 D up
                      %nbr_block1(2,4) will contain block1 D down
    nbr_surr1=zeros(1,4);

    nbr_block2=zeros(2,16);%nbr_block2(1,1) will contain block2 AA up
                      %nbr_block1(2,1) will contain block2 AA down
                      % . . . 
                      %nbr_block1(1,5) will contain block2 BA up
                      %nbr_block1(2,5) will contain block1 BA down
                      
    nbr_surr2=zeros(1,16);

    nbr_block3=zeros(2,64);
    nbr_surr3=zeros(1,64);


    nbr_block4=zeros(2,256);
    nbr_surr4=zeros(1,256);


    P_block1=zeros(2,4);%in row 1 prob of having block j and spin up
                    %in row 2 prob of having block j and spin down

    P_block2=zeros(2,16);

    P_block3=zeros(2,64);

    P_block4=zeros(2,256);



    P_given_b1=zeros(2,4);
    P_given_b2=zeros(2,16);
    P_given_b3=zeros(2,64);
    P_given_b4=zeros(2,256);


for a=1:grid_dim
    
    for i=1:4
        if spin_list(list_block_dim1(a,:))==block1Type(:,i)  
            nbr_surr1(i)=nbr_surr1(i)+1;
            if spin_list(a)==1
                nbr_block1(1,i)=nbr_block1(1,i)+1;
            else
                nbr_block1(2,i)=nbr_block1(2,i)+1;
            end
        end
        
    end
    
    for j=1:16
        if spin_list(list_block_dim2(a,:))==block2Type(:,j)  
            nbr_surr2(j)=nbr_surr2(j)+1;
            if spin_list(a)==1
                nbr_block2(1,j)=nbr_block2(1,j)+1;
            else
                nbr_block2(2,j)=nbr_block2(2,j)+1;
            end
        end
        
    end
    
    
    for k=1:64
        if spin_list(list_block_dim3(a,:))==block3Type(:,k)  
            nbr_surr3(k)=nbr_surr3(k)+1;
            if spin_list(a)==1
                nbr_block3(1,k)=nbr_block3(1,k)+1;
            else
                nbr_block3(2,k)=nbr_block3(2,k)+1;
            end
        end
        
    end
    
    
    for v=1:256
        
        if spin_list(list_block_dim4(a,:))==block4Type(:,v)  
            nbr_surr4(v)=nbr_surr4(v)+1;
            if spin_list(a)==1
                nbr_block4(1,v)=nbr_block4(1,v)+1;
            else
                nbr_block4(2,v)=nbr_block4(2,v)+1;
            end
        end
        
    end
    
     
    
end


  magnetization(stat)=sum(spin_list)/grid_dim;


  P_block1(1,:)=nbr_block1(1,:)./grid_dim;
  P_block1(2,:)=nbr_block1(2,:)./grid_dim;


  P_block2(1,:)=nbr_block2(1,:)./grid_dim;
  P_block2(2,:)=nbr_block2(2,:)./grid_dim;


  P_block3(1,:)=nbr_block3(1,:)./grid_dim;
  P_block3(2,:)=nbr_block3(2,:)./grid_dim;


  P_block4(1,:)=nbr_block4(1,:)./grid_dim;
  P_block4(2,:)=nbr_block4(2,:)./grid_dim;


  P_surr1=nbr_surr1./grid_dim;
  P_surr2=nbr_surr2./grid_dim;
  P_surr3=nbr_surr3./grid_dim;
  P_surr4=nbr_surr4./grid_dim;
  
  
%maybe create a function where you send as input the vector type of block
%and the state of the system and take as output the porbability of having
%x up and down.

%% COMPUTE BLOCK ENTROPY


for i=1:4
  if P_surr1(i)==0
     P_given_b1(1,i)=0;
     P_given_b1(2,i)=0;
  else
      P_given_b1(1,i)=P_block1(1,i)/P_surr1(i);
      P_given_b1(2,i)=P_block1(2,i)/P_surr1(i);
  end
end


for j=1:16
  if P_surr2(j)==0
     P_given_b2(1,j)=0;
     P_given_b2(2,j)=0;
  else
      P_given_b2(1,j)=P_block2(1,j)/P_surr2(j);
      P_given_b2(2,j)=P_block2(2,j)/P_surr2(j);
  end
end


for k=1:64
  if P_surr3(k)==0
     P_given_b3(1,k)=0;
     P_given_b3(2,k)=0;
  else
      P_given_b3(1,k)=P_block3(1,k)/P_surr3(k);
      P_given_b3(2,k)=P_block3(2,k)/P_surr3(k);
  end
end




for v=1:256
  if P_surr4(v)==0
     P_given_b4(1,v)=0;
     P_given_b4(2,v)=0;
  else
      P_given_b4(1,v)=P_block4(1,v)/P_surr4(v);
      P_given_b4(2,v)=P_block4(2,v)/P_surr4(v);
  end
end



S_block1(stat)=0;


for i=1:4
    for k=1:2
        if P_given_b1(k,i)~=0  
        S_block1(stat)=S_block1(stat)+P_surr1(i)*P_given_b1(k,i)*log(1/(P_given_b1(k,i)));
        end
    end
end

S_block2(stat)=0;


for j=1:16
    for k=1:2
        if P_given_b2(k,j)~=0  
        S_block2(stat)=S_block2(stat)+P_surr2(j)*P_given_b2(k,j)*log(1/(P_given_b2(k,j)));
        end
    end
end

S_block3(stat)=0;

for n=1:64
    for k=1:2
        if P_given_b3(k,n)~=0  
        S_block3(stat)=S_block3(stat)+P_surr3(n)*P_given_b3(k,n)*log(1/(P_given_b3(k,n)));
        end
    end
end
 
S_block4(stat)=0;

for v=1:256
    for k=1:2
        if P_given_b4(k,v)~=0  
        S_block4(stat)=S_block4(stat)+P_surr4(v)*P_given_b4(k,v)*log(1/(P_given_b4(k,v)));
        end
    end
 end




% iter
%   for i=1:grid_dim
%       if spin_list(i)==1
%           plot(grid_pos(i,1),grid_pos(i,2),'r*');hold on
%       else plot(grid_pos(i,1),grid_pos(i,2),'b*');hold on 
%       end
%   end
%   drawnow;

end%after burn in


end%end of total loop

S_block1_av(temp)=sum(S_block1)/stat;
S_block2_av(temp)=sum(S_block2)/stat;
S_block3_av(temp)=sum(S_block3)/stat;
S_block4_av(temp)=sum(S_block4)/stat;



magn_av=sum(magnetization)/stat;

if magn_av==1
    S_magn(temp)=0;
elseif magn_av==-1
    S_magn(temp)=0;
else 
S_magn(temp)=-((1+magn_av)/2*log((1+magn_av)/2)+(1-magn_av)/2*log((1-magn_av)/2));
end

%% SAVE CONFIGURATION 

conf_T(:,temp)=spin_list;
figure 

  for i=1:grid_dim
      if conf_T(i,temp)==1
           plot(grid_pos(i,1),grid_pos(i,2),'r*');hold on
       else plot(grid_pos(i,1),grid_pos(i,2),'b*');hold on 
       end
   end



end %end of temperature confront



figure
plot(T,S_block1_av,'k*',T,S_block2_av,'b*',T,S_block3_av,'c*',T,S_block4_av,'g*');
xlabel('T');
ylabel('S');
grid on;



