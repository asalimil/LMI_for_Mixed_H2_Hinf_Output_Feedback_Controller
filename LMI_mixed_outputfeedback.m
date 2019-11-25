clc; clear; close all;

A=[-1 1 0 1 0 1;-1 -2 -1 0 0 1;1 0 -2 -1 1 1;-1 1 -1 -2 0 0;-1 -1 1 1 -3 -1;0 -1 0 0 -1 -2];
B=[0 -1 -1;0 0 0;-1 1 1;-1 0 0;0 0 1;-1 1 1];
C=[0 1 0 -1 -1 -1;0 0 0 -1 0 0;1 0 0 0 -1 0];
D=[0 0 0;0 0 0;0 0 0];

% Extracting the 9-matrix
B1=[B zeros(6,3)];
B2=B;
C1=[C;zeros(3,6)];
C2=C;
D11=[D zeros(3);zeros(3) zeros(3)];
D12=[D;eye(3)];
D21=[D eye(3)];
D22=D;
P=[A B1 B2; C1 D11 D12; C2 D21 D22];

% Represent the P matrix in A B C D format to be used in sys command
At=A;
Bt=[B1 B2];
Ct=[C1;C2];
Dt=[D11 D12;D21 D22];
sys1=ss(At,Bt,Ct,Dt);

% LMIs
beta1=sdpvar(1);
beta2=sdpvar(1);
X1=sdpvar(6);
Y1=sdpvar(6);
Z=sdpvar(6);
An=sdpvar(6,6);
Bn=sdpvar(6,3,'full');
Cn=sdpvar(3,6,'full');
Dn=sdpvar(3,3);
Const=[];
M1=[A*Y1+Y1*A'+B2*Cn+Cn'*B2'      (A'+An+(B2*Dn*C2)')'               (B1+B2*Dn*D21); 
        (A'+An+(B2*Dn*C2)')           X1*A+A'*X1+Bn*C2+C2'*Bn'           (X1*B1+Bn*D21);
        (B1+B2*Dn*D21)'                           (X1*B1+Bn*D21)'                                     -eye(6)];
    
M2=[Y1                                       eye(6)                (C1*Y1+D12*Cn)';
        eye(6)                                     X1                  (C1+D12*Dn*C2)';
       (C1*Y1+D12*Cn)       (C1+D12*Dn*C2)                                Z];
   
Const=[Const;  (D11+D12*Dn*D21)==0; trace(Z) <= beta1];

M3=[A*Y1+Y1*A'+B2*Cn+Cn'*B2'   (A'+An+(B2*Dn*C2)')'   (B1+B2*Dn*D21)           (C1*Y1+D12*Cn)'; 
        (A'+An+(B2*Dn*C2)')   X1*A+A'*X1+Bn*C2+C2'*Bn'     (X1*B1+Bn*D21)         (C1+D12*Dn*C2)';
        (B1+B2*Dn*D21)'              (X1*B1+Bn*D21)'               -beta2*eye(6)                  (D11+D12*Dn*D21)';
         (C1*Y1+D12*Cn)             (C1+D12*Dn*C2)              (D11+D12*Dn*D21)                             -eye(6)];
     
Const=[Const; M1 <= 0];
Const=[Const; M2 >= 0];
Const=[Const; M3 <= 0];
beta=beta1+beta2;
optimize(Const,beta);
beta=value(beta);
value(beta1);
value(beta2);
gamma1=sqrt(value(beta1));
gamma2=sqrt(value(beta2));
H_2_optimal_gain=value(gamma1)
H_infinity_optimal_gain=value(gamma2)
