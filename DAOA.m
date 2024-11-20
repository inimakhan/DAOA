%__________________________________________________________________ %
%                                                                   %
%       Dynamic Arithmetic Optimization Algorithm  (DAOA)           %
%                                                                   %
%                                                                   %
%               Developed in MATLAB R2020b (MacOs-Big Sur)          %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%       Seyedali Mirjalili    (ʘ‿ʘ)     Nima Khodadadi              %
%                                                                   %
%                                                                   %
%                                                                   %
%                            e-Mail(1)                              %
%                ---------------------------------                  %
%                    ali.mirjalili@gmail.com                        %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%                            e-Mail(2)                              %
%                ---------------------------------                  %
%                         inimakhan@me.com                          %
%                                                                   %
%                                                                   %
%                                                                   %
%                           Homepage(1)                             %
%                ---------------------------------                  %
%                    http://www.alimirjalili.com                    %
%                                                                   %
%                           Homepage(2)                             %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function Details

%  function [Best_FF,Best_P,Conv_curve]=DAOA(N,M_Iter,LB ,UB,Dim,ObjFuncName)
display('AOA Working');

Dim=2;                          % Number of Variable
LB = (0)*ones(1,Dim);        % Upper Bound
UB = (1)*ones(1,Dim);           % Lower Bound
F_obj= @F1;                    %Name of the Function
N=5;                           %Number of Colors (Npop)
M_Iter=500; 
Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);
%% Initialize the positions of solution

X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));          % (fitness values)
Ffun_new=zeros(1,size(Xnew,1));   % (fitness values)


C_Iter=1;
Mu=0.001;
alpha=25;

for i=1:size(X,1)   
    Ffun(1,i)=F_obj(X(i,:));       %Calculate the fitness values of solutions
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end

%% Main Loop

while C_Iter<M_Iter+1                        % Main loop
    
   %DAF=((C_Iter/M_Iter+1)^((-1)*alpha));    %DAF1
    DAF=(M_Iter+1/C_Iter)^((alpha));        %DAF2
    
    DCS=.99*((1-(C_Iter/M_Iter)^(0.5)));     % DCS
    
    
    %Upate the Position of solutions
    
    for i=1:size(X,1)                        % if each of the UB and LB has a just value
        for j=1:size(X,2)
            r1=rand();
            if (size(LB,2)==1)
                if r1<DAF
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=(Best_P(1,j)/(DCS+eps)*((UB-LB)*Mu+LB));
                    else
                        Xnew(i,j)=(Best_P(1,j)*DCS*((UB-LB)*Mu+LB));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=(Best_P(1,j)-DCS*((UB-LB)*Mu+LB));
                    else
                        Xnew(i,j)=(Best_P(1,j)+DCS*((UB-LB)*Mu+LB));
                    end
                end
            end
            
            
            if (size(LB,2)~=1)                          % if each of the UB and LB has more than one value
                r1=rand();
                if r1<DAF
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=((Best_P(1,j)/(DCS+eps)*((UB(j)-LB(j))*Mu+LB(j))));
                    else
                        Xnew(i,j)=((Best_P(1,j)*DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=((Best_P(1,j)-DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    else
                        Xnew(i,j)=((Best_P(1,j)+DCS*((UB(j)-LB(j))*Mu+LB(j))));
                    end
                end
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB;                    % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB;                    % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
        
        Ffun_new(1,i)=F_obj(Xnew(i,:));          % calculate Fitness function
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
            Best_FF=Ffun(1,i);
            Best_P=X(i,:);
        end
        
    end
    
    
    %Update the convergence curve
    
    Conv_curve(C_Iter)=Best_FF;
    plot(Conv_curve,'Color','r','LineWidth',2)
    title('Convergence curve')
    xlabel('Iteration');
    ylabel('Best fitness function');
    axis tight
    legend('DAOA')
    %Print the best solution details after every 50 iterations
    %     if mod(C_Iter,50)==0
    display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
    %     end
    
    C_Iter=C_Iter+1;  % incremental iteration
    
 end

