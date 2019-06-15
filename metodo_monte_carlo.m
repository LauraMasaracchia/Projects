%ciclo for che mi calcola i risultati un po' di volte
tot= input('quanti confronti in temperatura? ')

N=input('scegli la lunghezza del reticolo N = ')
n= input('numero interazioni che vuoi fare: n = ')
J=input('stabilisci coefficiente di interazione J = ')
H=input('campo magnetico esterno (in tesla) H = ')

Kb=1.380;

temp=zeros(1,tot);
magnet=zeros(1,tot);

for w=1:tot
    w
    temp(w)=input('temperatura= ')
    
    B=1/(Kb*temp(w));
    
    MATRIX=zeros(N,N);
   s=zeros(1,4);
%dovresti avere una matrice riempita casualmente di meno uno e uno
   for i=1:N
      for j=1:N
         r=rand;
           if r<=0.5
              MATRIX(i,j)=1;
           end
           if r>0.5
              MATRIX(i,j)=-1;
           end
      end
   end
    
  %scegliamo a caso un elemento della matrice: come fare? magari 
  %scegliendo un numero random da convertire in un elemento della matrice

 %r sarà un num compreso tra 0 e 1 
 %prendiamo a caso coordinate x e y
   for j=1:n
    
      p1=rand*N;
     for i=1:N
         if i>=p1 && i<p1+1
            p1=i;
          break
         end
     end
       p2=rand*N;
      for i=1:N
          if i>=p2 && i<p2+1
             p2=i;
           break
          end
      end
     p1; 
     p2;

     x=MATRIX(p1,p2); %il nostro punto
       if p1==1 
          p1u=N;
       else p1u=p1-1;
       end
      if p2==1
         p2u=N;
      else p2u=p2-1;
      end
       
      if p1==N
         p1d=1;
      else p1d=p1+1;
      end    
      if p2==N
         p2d=1;
      else p2d=p2+1;
      end
   
     s(1)=MATRIX(p1u,p2);
     s(2)=MATRIX(p1,p2u);
     s(3)=MATRIX(p1d,p2);
     s(4)=MATRIX(p1,p2d);
     S=sum(s,2);
      E1= -J*S - H; %energia del vicinato se il nostro m(j) diventa -1
      E2=J*S + H; %se il nostro m(j) diventa -1
      fun1=(exp(-B*E1))/((exp(-B*E1))+(exp(-B*E2)));
      r=rand;
      j;
     if r<=fun1
         MATRIX(p1,p2)=1;
     else
        MATRIX(p1,p2)=-1;
     end
   MATRIX;
   
   end

  prov=sum(MATRIX);
  M=sum(prov);
  m=M/(N^2)
  magnet(w)=m;


   fup=zeros(N,N);
   fdown=zeros(N,N);
   
   for t=1:N
       for r=1:N
           if MATRIX(t,r)==1
               fup(t,r)=r;
           else fdown(t,r)=r;
           end
       end
   end
plot (fup,'ro');
hold on
plot(fdown,'bo');
pause
close all



end



   z=linspace(0,max(temp));
  % c=polyfit(temp,magnet,tot-1);
  % v=polyval(c,z);
plot(temp,magnet,'og');
%plot(temp, magnet,'r');

