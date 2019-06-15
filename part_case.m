close all
clear all

% INITIALIZATION of the system. 

%setting J to 1 aand Kb to 1 the transition temperature is around 2 (?)

Kb=1;
T=[0.5];
J=1;


burn_in=2000000;%iterations to reach equilibrium, empirically set
indep_stat=50000;
total_nbr_iter=burn_in+(indep_stat*20);
beta=1./(Kb*T);
side_dim=79;
grid_dim=side_dim*side_dim;

grid_pos=zeros(grid_dim,2);%list of positions of each grid site
spin_list=zeros(grid_dim,1);%list of spin of each grid site
conf_T=zeros(grid_dim,length(T));


 a=1;

for i=1:side_dim
    for j=1:side_dim
        grid_pos(a,:)=[j,i];
        a=a+1;
    end
end

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



%block neighbors

list_block_dim1=zeros(grid_dim,2);
list_block_dim2=zeros(grid_dim,4);
list_block_dim3=zeros(grid_dim,6);
list_block_dim4=zeros(grid_dim,8);

%list containing in each row a the numbers of the position forming the block
%of dimension one around position a



for a=1:grid_dim
    nl=nn_list(a,2);
    nd=nn_list(a,4);
    nll=nn_list(nn_list(a,2),2);
    ndr=nn_list(nn_list(a,4),1);
    nlll=nn_list(nn_list(nn_list(a,2),2),2);
    ndrr=nn_list(nn_list(nn_list(a,4),1),1);
    nllll=nn_list(nn_list(nn_list(nn_list(a,2),2),2),2);
    ndrrr=nn_list(nn_list(nn_list(nn_list(a,4),1),1),1);
    
list_block_dim1(a,:)=[nl,nd]; %left, down
list_block_dim2(a,:)=[nl,nd,nll,ndr]; %left, down, leftleft, downright
list_block_dim3(a,:)=[nl,nd,nll,ndr,nlll,ndrr]; %left, down, leftleft, downright, leftleftleft, downrightright
list_block_dim4(a,:)=[nl,nd,nll,ndr,nlll,ndrr,nllll,ndrrr];
end


%for any element x indicated in row a we will consider a block1 of the kind
%
%              
% o3  o2  o1  x
%             o1  o2  o3
%






block1Type=zeros(2,4);%in each column a kind of block. need it to be in a
%column since it has to be compared with the spin list, which is a column.

%left
%down

block1Type(:,1)=[1;1];%A
block1Type(:,2)=[-1;-1];%B
block1Type(:,3)=[-1;1];%C
block1Type(:,4)=[1;-1];%D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example: block1,type A
%     A           AA           AD           AB           AC
%                                                      
%  1  x       1 1 x       1 1 x        -1 1 x        -1 1 x  
%     1           1 1         1 -1          1 -1          1 1
%
%for all the x count how many blocks of type 1 are there and when x=1 or -1
%Then make a statistic.



block2Type=zeros(4,16);
%left
%down
%leftleft
%downright
c=1;
for i=1:4
    for j=1:4
        block2Type(:,c)=[block1Type(:,i);block1Type(:,j)];%AA,AB,AC,AD,BA,BB,BC,BD,CA,CB,CC,CD,DA,DB,DC,DD
            c=c+1;
    end
end


block3Type=zeros(6,64);
%left
%down
%leftleft
%downright
%leftleftleft
%downrightright

c=1;
for i=1:16
    for j=1:4
        block3Type(:,c)=[block2Type(:,i);block1Type(:,j)];%AAA,AAB,...,AAD,ABA,...ABD,ACA...,ACD,...,DDA,...,DDD
            c=c+1;
    end
end


block4Type=zeros(8,256);
%left
%down
%leftleft
%downright
%leftleftleft
%downrightright
%leftleftleftleft
%downrightrightright

c=1;
for i=1:64
    for j=1:4
        block4Type(:,c)=[block3Type(:,i);block1Type(:,j)];
            c=c+1;
    end
end




%% pattern configuration


% lineS  1 0 1 0 1
%        1 0 1 0 1
%        1 0 1 0 1
%
%block entropy=0, magnetization entropy= log2

%NOTE THAT WHEN WE HAVE ODD NUMBER OF ATOMS IN A ROW THIS BECOMES
%  1 0 1 0 1
%  0 1 0 1 0  
%  1 0 1 0 1
%
% AND IT STILL SPOTS SOME RESIDUAL ENTROPY WITH THE BLOCKS
% magnetization entropy= log2


 for a=1:grid_dim
     if mod(a,2)==0 %is even
     spin_list(a)=+1;
     else
         spin_list(a)=-1;
     end
 end 



%new pattern: lines of bigger blocks
%   1 1 1 0 0 0 1 1 1
%   1 1 1 0 0 0 1 1 1
%   1 1 1 0 0 0 1 1 1
%blocks entropy zero
%
% in case of non even number of atoms per row the pattern is not a stripe
% but diagonals. Here the block entropy is something but goes down as soon
% as we encrease the block dimension


% 
% tt=0;
% for a=1:grid_dim
%     if tt==0 || tt==1 || tt==2
%         spin_list(a)=1;
%         tt=tt+1;
%     else spin_list(a)=-1;
%         tt=tt+1;
%     end
%     if tt==6
%         tt=0;
%     end
% end



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



S_block1=0;


for i=1:4
    for k=1:2
        if P_given_b1(k,i)~=0  
        S_block1=S_block1+P_surr1(i)*P_given_b1(k,i)*log(1/(P_given_b1(k,i)));
        end
    end
end

S_block2=0;


for j=1:16
    for k=1:2
        if P_given_b2(k,j)~=0  
        S_block2=S_block2+P_surr2(j)*P_given_b2(k,j)*log(1/(P_given_b2(k,j)));
        end
    end
end

S_block3=0;

for n=1:64
    for k=1:2
        if P_given_b3(k,n)~=0  
        S_block3=S_block3+P_surr3(n)*P_given_b3(k,n)*log(1/(P_given_b3(k,n)));
        end
    end
end
 
S_block4=0;

for v=1:256
    for k=1:2
        if P_given_b4(k,v)~=0  
        S_block4=S_block4+P_surr4(v)*P_given_b4(k,v)*log(1/(P_given_b4(k,v)));
        end
    end
end




  for i=1:grid_dim
       if spin_list(i)==1
           plot(grid_pos(i,1),grid_pos(i,2),'r*');hold on
       else plot(grid_pos(i,1),grid_pos(i,2),'b*');hold on 
       end
  end

magn_av=sum(spin_list)/grid_dim;
   
if magn_av==1
    S_magn=0;
elseif magn_av==-1
    S_magn=0;
else 
S_magn=-((1+magn_av)/2*log((1+magn_av)/2)+(1-magn_av)/2*log((1-magn_av)/2));
end

  
  









    
    
    
    

