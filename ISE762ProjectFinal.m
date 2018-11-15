%% Common Random Numbers 
clear all; 
clc; 
clf; 
tic; 
lambda = 1/6; 
mu1 = 1/15;
muk = 1/20;
mu2 = 1/5; 
p = .75;

N = 500; 
NSim = 200; % Number of simulation runs

NS = 3; %Number of servers at station 1
NK = 5; %Number of kiosks


%Uncomment if we want to vary parameters
%for NK = 3:20
%for lambda = .1:.01:.2
    
W1 = zeros(N,1); % Initialize waiting times vector for the queues
W2 = zeros(N,1);

PTk = [];
PTs = [];
PTQ1 = [];
SJ1 = []; %Sojourn time for Alt 1
SJ2 = []; %Sojourn time for Alt 2

for k=1:NSim 
    
    %Generate interarrival times and service times 
    A = zeros(N,1); %Arrival times
    S1 = zeros(N,1); %Service times at station 1
    SK = zeros(N,1); %Service times at kiosks
    S2 = zeros(N,1); %Service times for station 2

    %Initialize values
    A(1) = -log(rand)/lambda;
    A(N+1) = Inf; 
    S1(1) = -log(rand)/mu1;
    SK(1) = -log(rand)/muk;
    S2(1) = -log(rand)/mu2;
    for i=2:N
        A(i) = A(i-1)-log(rand)/lambda;
        S1(i) = -log(rand)/mu1;
        SK(i) = -log(rand)/muk; 
        S2(i) = -log(rand)/mu2;
    end 
    
    t = 0; 
    NA = 0;
    NDS = 0; 
    ND2 = 0; 
    S = Inf(NS,1); %Departure times from servers
     
    tD2 = Inf;
    
    Output1 = []; 
    tA = A(1); %Set time of first arrival 
    QS = []; 
    Q1 = [];
    Q2 = []; 
    SID = 1:NS;

    
    Bs = [];
    while ND2 < N
        
        if tA <= min([S;tD2])
            t = tA; 
            NA = NA + 1;
            if sum(S<Inf) < NS %if number of customers in system is less than the number of servers
                idx = find(S==Inf,1);
                QS(idx) = NA; 
                S(idx) = t + S1(NA);
                Bs = [Bs; NA SID(idx) t S(idx)];
            else %all servers full 
                QS(end+1) = NA; 
            end 
            if NA < N
                tA = A(NA+1);
            else 
                tA = Inf;
            end 
            Output1 = [Output1; NA t 0 0]; 
            
        elseif min(S) <= tD2
            [t, I] = min(S);
            Output1(QS(I),3) = t;
            NDS = NDS + 1;
            
            
            Q2 = [Q2 QS(I)];
            if length(Q2) == 1
                tD2 = t + S2(Q2(1));
            end
            
            if length(QS) > NS %If there were more people than kiosks
                QS(I) = QS(NS+1); %Next person in line goes to the vacant kiosk
                QS(NS+1) = []; %Vacate the spot they left
                S(I) = t + S1(QS(I)); %Assign departure time
                
                Bs = [Bs; QS(I) SID(I) t S(I)];
            else
                QS(I) = []; 
                S(I) = [];
                S(end+1) = Inf; %Mark kiosk as vacant
                
                temp = SID(I);
                SID(I) = [];
                SID(end+1) = temp;
            end 
        else 
            t = tD2;
            ND2 = ND2 + 1;
            Output1(Q2(1),4) = t;
            if length(Q2) == 1 
                tD2 = Inf; 
                Q2 = [];
            else
                tD2 = t + S2(Q2(1));
                Q2 = Q2(2:end);
            end
        end 

        
    end 
    
    %Compute the waiting times
    pts = zeros(1,NS);
    for a = 1:NS
        idxx  = Bs(Bs(:,2)==a);
        pts(a) = sum(Bs(idxx,4) - Bs(idxx,3))/max(Output1(:,3));
    end
    PTs = [PTs; pts];

    w = Output1(:,4) - Output1(:,2); % total waiting time
    W1 = W1 + w;

    SJ1 = [SJ1 max(Output1(:,4))-Output1(1,2)];


%% System 2 





    t = 0; 
    NA = 0;
    NDK = 0; 
    ND1 = 0; 
    ND2 = 0; 
    K = Inf(1,NK); %Departure times from kiosks
    tD1 = Inf; 
    tD2 = Inf;
    
    Output2 = []; 
    tA = A(1); %Set time of first arrival 
    QK = []; 
    Q1 = [];
    Q2 = [];
    BQ1 = zeros(1,N);
    
    Bk = [];
    KID = 1:NK;
    KID = KID';
    
    while ND2 < N
        
        if tA <= min([K tD1 tD2])
            t = tA; 
            NA = NA + 1;
            if sum(K<Inf) < NK %if number of customers in system is less than the number of kiosks
                idx = find(K==Inf,1);
                QK(idx) = NA; 
                K(idx) = t + SK(NA);
                
                Bk = [Bk; NA KID(idx) t K(idx)];
            else %all kiosks full 
                QK(end+1) = NA; 
            end 
            if NA < N
                tA = A(NA+1);
                
            else 
                tA = Inf;
            end 
            Output2 = [Output2; NA t 0 0 0]; 
            
        elseif min(K) <= min(tD1,tD2)
            [t, I] = min(K);
            Output2(QK(I),3) = t;
            NDK = NDK + 1;
            
            if rand <= p %If kiosk was able to handle
                Q2 = [Q2 QK(I)]; %Go to queue 2
                if length(Q2) == 1
                    tD2 = t + S2(Q2(1));
                end
            else %Else, go to queue 1
                Q1 = [Q1 QK(I)];
                BQ1(QK(I)) = 1; %Mark that this customer had to use Q1
                if length(Q1) == 1 %If Q1 was empty
                    tD1 = t + S1(Q1(1)); %assign departure time
                end 
            end
            if length(QK) > NK %If there were more people than kiosks
                QK(I) = QK(NK+1); %Next person in line goes to the vacant kiosk
                QK(NK+1) = []; %Vacate the spot they left
                K(I) = t + SK(QK(I)); %Assign departure time
                
                Bk = [Bk; QK(I) KID(I) t K(I)];
            else
                QK(I) = []; 
                K(I) = [];
                K(end+1) = Inf; %Mark kiosk as vacant
                
                temp = KID(I);
                KID(I) = [];
                KID(end+1) = temp;
            end 
        elseif tD1 <= tD2
            t = tD1;
            ND1 = ND1 + 1;
            Output2(Q1(1),4) = t;
            Q2 = [Q2 Q1(1)];
            if length(Q1) == 1 %If the queue had only 1 person
                tD1 = Inf; 
                Q1 = []; %Empty the queue
            else
                tD1 = t + S1(Q1(1));
                Q1 = Q1(2:end);
            end
            if length(Q2) == 1
                tD2 = t + S2(Q2(1)); 
            end 
            
        else 
            t = tD2;
            ND2 = ND2 + 1;
            Output2(Q2(1),5) = t;
            if length(Q2) == 1 
                tD2 = Inf; 
                Q2 = [];
            else
                tD2 = t + S2(Q2(1));
                Q2 = Q2(2:end);
            end
        end 
        
    end 
    
    %Compute the waiting times
    
    w = Output2(:,5) - Output2(:,2); % total waiting time
    W2 = W2 + w;
    ptk = zeros(1,NK);
    ptb = 0;
    for a = 1:NK
        idxx  = Bk(Bk(:,2)==a);
        ptk(a) = sum(Bk(idxx,4) - Bk(idxx,3))/max(Output2(:,3));
    end
    PTk = [PTk; ptk];
    ptq1 = sum(BQ1.*S1')/max(Output2(:,4));
    PTQ1 = [PTQ1; ptq1]; 
    SJ2 = [SJ2 max(Output2(:,5))-Output2(1,2)];
end


EW1 = W1/NSim;
EW2 = W2/NSim;

figure 
plot((1:1:N),EW1-EW2)
title(('Reduction in Waiting Time from Original System to Kiosk Check-In System'))
xlabel('Customer')
ylabel('Time Saved (minutes)')