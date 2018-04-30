 //Projet d'optimisation 
//Kenneth et Anne-Marie

//Question 1

//Utilisation de fscanfMat
N=fscanfMat("C:\Users\anne\Documents\"+"profil.txt");
//Représentation de la courbe (xi,gi)i
X0=N(:,1);
G0=N(:,2);
clf
plot(X0,G0,'b')
xlabel("Points de la route","frontsize",6,"color","green")
ylabel("Altitude", "fontsize", 6, "color", "blue");

//Question 2

 //Récupération dqns X0 et G0
n=185
for i=1:n
   X(i)=X0(i);  
    end
 for i=1:n
    G(i)=G0(i);
 end

//Calcul de la pente locale maximale
h=(X(2)-X(1));
p=0;
for i=1:184
    tmp=(G0(i+1)-G0(i))/h;
   if tmp>p
       p=tmp;
    end
     
end
printf("La pente locale maximale est %f\n",p);
alpha=0.1;
//Question 3

//Definition des matrices A et C et du vec

v=ones(1,n);
v(1)=0.5;
v(n)=0.5;
A=h*diag(v);

C1=zeros(n-1,n);
for i=1:n-1
    for j=1:n
        if i==j then
            C1(i,j)=-1;
            C1(i,j+1)=1
         else
         end
     end
        
end

b=alpha*h*ones(2*(n-1),1);

C2=zeros(n-1,n);
for i=1:n-1
    for j=1:n
        if i==j then
            C2(i,j)=1;
            C2(i,j+1)=-1
         else
         end
     end
        
end
C=[C1;C2];
//Question 4 inverse de la matrice A
w=zeros(1,n);
w(1)=2/h;
w(n)=2/h;
for i=2:(n-1)
    w(i)=1/h;
end
Am1=diag(w);
//Question 5
//Programmation de l'Algorithme d'Uzawa

function [U1,z,k]=minimize_uzawa(A,G,C,b,U0,z0,rh0,eps,kmax)
    k=1;
    U1=0.5*Am1*(2*A*G-C'*z0);
    z=max (z0+rh0*(C*U1-b),0);
    r=2*A*G; 
    while norm(U1-U0)>=eps & k<kmax
                U0=U1;
                U1=0.5*Am1*(r-C'*z);
                z=max(z+rh0*(C*U1-b),0);
                k=k+1;
             
          plot(X,U1,'r',X,G);
 
    end
endfunction
funcprot(0)
// Programme principale
eps=0.01;
kmax=20000;
U0=2*ones(n,1);
z0=ones(2*(n-1),1)
rh0=0.33;
[U1,z,k]=minimize_uzawa(A,G,C,b,U0,z0,rh0,eps,kmax)
         
//Question 6        
         
 figure 0
  clf     
  plot(X,U1,'r',X,G);
  legend(['Tracé final du terrain';'Tracé  initial du terrain'])
 xtitle("Points de discrétisation du terrain")
           
printf('Uzawa:avec %i itération(s)\n',k);

//Question 7   Variation de alpha
//Se fera depuis le haut
//On remarque que quand alpa diminue,le nombre d'tération augmente et vis-versa

//Question 8 observation si alpa>p on a un petit nombre d'itérations

//Question 9 
alph=[0.33,0.22,0.1,0.01];
A
C
G
eps=0.01;
kmax=20000;
U0=2*ones(n,1);
z0=ones(2*(n-1),1)
rh0=0.33;
clf
for j=1:4
    b=alph(j)*h*ones(2*(n-1),1);
    [U1,z,k]=minimize_uzawa(A,G,C,b,U0,z0,rh0,eps,kmax)
    plot(X,U1,'g')
    end
 
//L'entrepreneur verra la pente alpha pour laquelle la courbe est la plus lisse

//Question 10
// Sur la solution,d'après la figure 1 qu'on a retiré plus de terrain qu'on en a remis

//Quesxtion 11
//Dès qu'on ajoute une contrainte supplé   //mentaire permettant de limiter ces angles//la méthode converge vers la solution plus//rapidement

//Question 12
