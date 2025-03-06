clear all;
close all;
clc;

M1 = [-3 -4]';
M2 = [9,8]';
M3 = [7,-3]';
M4 = [0, 8]';
S1 = [1.2 0.5; 0.5 1.2];
S2 = [0.6 0.2; 0.2 0.6];
S3 = [1 -0.7; -0.7 1];
S4 = [1.3 0.5; 0.5 1.3];

N = 500;K = 4;
X1 = mvnrnd(M1,S1,N)';
X2 = mvnrnd(M2,S2,N)';
X3 = mvnrnd(M3,S3,N)';
X4 = mvnrnd(M4,S4,N)';

figure();
plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bo', X3(1,:), X3(2,:), 'go', X4(1,:), X4(2,:), 'yo');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4');
xlabel('x1');ylabel('x2');
grid on;

%% C-mean


ll = zeros(1,10);
for k =1:10
pom = rand(1,4*N);
X11 = X1(:,1:20);
X21 = X2(:,1:20);
X31 = X3(:,1:20);
X41 = X4(:,1:20);


for i=20:N
    if pom(i)<0.25
        X11 = [X11 X1(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X21 = [X21 X1(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X31 = [X31 X1(:,i)];
    else
        X41 = [X41 X1(:,i)];
    end
end
for i=20:N
    if pom(i)<0.25
        X11 = [X11 X2(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X21 = [X21 X2(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X31 = [X31 X2(:,i)];
    else
        X41 = [X41 X2(:,i)];
    end
end
for i=20:N
    if pom(i)<0.25
        X11 = [X11 X3(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X21 = [X21 X3(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X31 = [X31 X3(:,i)];
    else
        X41 = [X41 X3(:,i)];
    end
end
for i=20:N
    if pom(i)<0.25
        X11 = [X11 X4(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X21 = [X21 X4(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X31 = [X31 X4(:,i)];
    else
        X41 = [X41 X4(:,i)];
    end
end

%{
figure();
plot(X11(1,:), X11(2,:), 'ro', X21(1,:), X21(2,:), 'bo', X31(1,:), X31(2,:), 'go', X41(1,:), X41(2,:), 'yo');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4');
xlabel('x1');ylabel('x2');
grid on;
title('Početna klasterizacija');
%}

N1 = length(X11);
N2 = length(X21);
N3 = length(X31);
N4 = length(X41);

M1 = mean(X11,2);
M2 = mean(X21,2);
M3 = mean(X31,2);
M4 = mean(X41,2);

lmax = 20;
reklas = 1;
l = 1;


while (l<lmax) && reklas
    X1pom = [];
    X2pom = [];
    X3pom = [];
    X4pom = [];
    reklas = 0;

    for i=1:N1
        d1 = sum((X11(:,i)-M1).^2);
        d2 = sum((X11(:,i)-M2).^2);
        d3 = sum((X11(:,i)-M3).^2);
        d4 = sum((X11(:,i)-M4).^2);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X11(:,i)];
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X11(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X11(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X11(:,i)];
            reklas = 1;
        end
    end
    for i=1:N2
        d1 = sum((X21(:,i)-M1).^2);
        d2 = sum((X21(:,i)-M2).^2);
        d3 = sum((X21(:,i)-M3).^2);
        d4 = sum((X21(:,i)-M4).^2);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X21(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X21(:,i)];
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X21(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X21(:,i)];
            reklas = 1;
        end
    end
    for i=1:N3
        d1 = sum((X31(:,i)-M1).^2);
        d2 = sum((X31(:,i)-M2).^2);
        d3 = sum((X31(:,i)-M3).^2);
        d4 = sum((X31(:,i)-M4).^2);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X31(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X31(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X31(:,i)];
        else
            X4pom = [X4pom X31(:,i)];
            reklas = 1;
        end
    end
    for i=1:N4
        d1 = sum((X41(:,i)-M1).^2);
        d2 = sum((X41(:,i)-M2).^2);
        d3 = sum((X41(:,i)-M3).^2);
        d4 = sum((X41(:,i)-M4).^2);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X41(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X41(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X41(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X41(:,i)];
        end
    end
    clear X11 X21 X31 X41
    X11 = X1pom;
    X21 = X2pom;
    X31 = X3pom;
    X41 = X4pom;
    
    %{
    figure();
    plot(X11(1,:), X11(2,:), 'ro', X21(1,:), X21(2,:), 'bo', X31(1,:), X31(2,:), 'go', X41(1,:), X41(2,:), 'yo');
    legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4');
    xlabel('x1');ylabel('x2');
    grid on;
    title(['Klasterizacija broj' num2str(l)]);
    pause
    %}

    N1 = length(X11);
    N2 = length(X21);
    N3 = length(X31);
    N4 = length(X41);

    M1 = mean(X11,2);
    M2 = mean(X21,2);
    M3 = mean(X31,2);
    M4 = mean(X41,2);
    l = l + 1;
end
ll(k) = l;
end
l_avg = mean(ll);

%% C-mean sa 2 predefinisana klastera
%i ovde i gore ako ne pripomognem na pocetku gore zabode a ovde iscrpi sve
%iteracije jer ima poneki odbirak koj uporno nece lepo da se klasifikuje pa
%im dajem pripomoc
ll = zeros(1,10);
for k =1:10
pom = rand(1,4*N);
X11 = [X1(:,1:10) X3(:,1:10)];
X21 = [X2(:,1:10) X4(:,1:10)];


for i=10:N
    if pom(i)<0.5
        X11 = [X11 X1(:,i)];
    else
        X21 = [X21 X1(:,i)];
    end
end
for i=10:N
    if pom(i)<0.5
        X11 = [X11 X2(:,i)];
    else
        X21 = [X21 X2(:,i)];
    end

end
for i=10:N
    if pom(i)<0.5
        X11 = [X11 X3(:,i)];
    else
        X21 = [X21 X3(:,i)];
    end
end
for i=10:N
    if pom(i)<0.5
        X11 = [X11 X4(:,i)];
    else
        X21 = [X21 X4(:,i)];
    end
end

%{

figure();
plot(X11(1,:), X11(2,:), 'ro', X21(1,:), X21(2,:), 'bo');
legend('Klasa 1', 'Klasa 2');
xlabel('x1');ylabel('x2');
grid on;
title('Početna klasterizacija');

%}

N1 = length(X11);
N2 = length(X21);

M1 = mean(X11,2);
M2 = mean(X21,2);

lmax = 20;
reklas = 1;
l = 1;


while (l<lmax) && reklas
    X1pom = [];
    X2pom = [];
    reklas = 0;

    for i=1:N1
        d1 = sum((X11(:,i)-M1).^2);
        d2 = sum((X11(:,i)-M2).^2);
        if d1<d2
            X1pom = [X1pom X11(:,i)];
        else
            X2pom = [X2pom X11(:,i)];
            reklas = 1;
        end
    end
    for i=1:N2
        d1 = sum((X21(:,i)-M1).^2);
        d2 = sum((X21(:,i)-M2).^2);
        if d1<d2
            X1pom = [X1pom X21(:,i)];
            reklas = 1;
        else
            X2pom = [X2pom X21(:,i)];
        end
    end
    clear X11 X21
    X11 = X1pom;
    X21 = X2pom;
    
    
    %{

    figure();
    plot(X11(1,:), X11(2,:), 'ro', X21(1,:), X21(2,:) ,'bo');
    legend('Klasa 1', 'Klasa 2');
    xlabel('x1');ylabel('x2');
    grid on;
    title(['Klasterizacija broj' num2str(l)]);
    pause
   
    %}
    

    N1 = length(X11);
    N2 = length(X21);

    M1 = mean(X11,2);
    M2 = mean(X21,2);
    l = l + 1;
end
ll(k) = l;
end
l_avg = mean(ll);



%% Metod kvadratne dekompozicije
close all;
clear all;
clc;

N = 500;
C = [6; 6];
R = 4*rand(1,N);
Teta = 2*pi*rand(1,N);
X = [R.*cos(Teta); R.*sin(Teta)];
X = X+C*ones(1,N);

C = [6; 8];
R = 8 + 2*rand(1,N);
Teta = pi*rand(1,N);
Y = [R.*cos(Teta); R.*sin(Teta)];
Y = Y+C*ones(1,N);

figure();
plot(X(1,:), X(2,:), 'ro', Y(1,:), Y(2,:) ,'bo');
legend('Klasa 1', 'Klasa 2');
xlabel('x1');ylabel('x2');
grid on;
title('Početne klase');

%% dve klase
ll = zeros(1,10);
for k =1:10
pom = rand(1,2*N);
X1 = [X(:,1:100)];
X2 = [Y(:,1:100)];


for i=100:N
    if pom(i)<0.5
        X1 = [X1 X(:,i)];
    else
        X2 = [X2 X(:,i)];
    end
end
for i=100:N
    if pom(i)<0.5
        X1 = [X1 Y(:,i)];
    else
        X2 = [X2 Y(:,i)];
    end

end

%{

figure();
plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bo');
legend('Klasa 1', 'Klasa 2');
xlabel('x1');ylabel('x2');
grid on;
title('Početna klasterizacija');

%}

N1 = length(X1);
N2 = length(X2);
P1 = N1/(2*N);
P2 = N2/(2*N);
M1 = mean(X1,2);
M2 = mean(X2,2);
S1 = cov(X1');
S2 = cov(X2');

lmax = 20;
reklas = 1;
l = 1;


while (l<lmax) && reklas
    X1pom = [];
    X2pom = [];
    reklas = 0;

    for i=1:N1
        d1 = 0.5*(X1(:,i)-M1)'*S1^(-1)*(X1(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X1(:,i)-M2)'*S2^(-1)*(X1(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        if d1<d2
            X1pom = [X1pom X1(:,i)];
        else
            X2pom = [X2pom X1(:,i)];
            reklas = 1;
        end
    end
    for i=1:N2
        d1 = 0.5*(X2(:,i)-M1)'*S1^(-1)*(X2(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X2(:,i)-M2)'*S2^(-1)*(X2(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        if d1<d2
            X1pom = [X1pom X2(:,i)];
            reklas = 1;
        else
            X2pom = [X2pom X2(:,i)];
        end
    end
    clear X1 X2
    X1 = X1pom;
    X2 = X2pom;
    
    
    %{

    figure();
    plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:) ,'bo');
    legend('Klasa 1', 'Klasa 2');
    xlabel('x1');ylabel('x2');
    grid on;
    title(['Klasterizacija broj' num2str(l)]);
    pause

    %}
    

    N1 = length(X1);
    N2 = length(X2);
    P1 = N1/(2*N);
    P2 = N2/(2*N);
    M1 = mean(X1,2);
    M2 = mean(X2,2);
    S1 = cov(X1');
    S2 = cov(X2');
    l = l + 1;
end
ll(k) = l;
end
l_avg = mean(ll);

%% cetiri klase

ll = zeros(1,10);
for k =1:10
pom = rand(1,4*N);
X1 = X(:,1:100);
X2 = X(:,101:200);
X3 = Y(:,1:100);
X4 = Y(:,101:200);


for i=200:N
    if pom(i)<0.25
        X1 = [X1 X(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X2 = [X2 X(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X3 = [X3 X(:,i)];
    else
        X4 = [X4 X(:,i)];
    end
end
for i=200:N
    if pom(i)<0.25
        X1 = [X1 Y(:,i)];
    elseif pom(i)>=0.25 && pom(i)<0.5
        X2 = [X2 Y(:,i)];
    elseif pom(i)>=0.5 && pom(i)<0.75
        X3 = [X3 Y(:,i)];
    else
        X4 = [X4 Y(:,i)];
    end
end


%{

figure();
plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bo', X3(1,:), X3(2,:), 'go', X4(1,:), X4(2,:), 'yo');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4');
xlabel('x1');ylabel('x2');
grid on;
title('Početna klasterizacija');

%}

N1 = length(X1);
N2 = length(X2);
N3 = length(X3);
N4 = length(X4);
P1 = N1/(2*N);
P2 = N2/(2*N);
P3 = N3/(2*N);
P4 = N4/(2*N);
M1 = mean(X1,2);
M2 = mean(X2,2);
M3 = mean(X3,2);
M4 = mean(X4,2);
S1 = cov(X1');
S2 = cov(X2');
S3 = cov(X3');
S4 = cov(X4');

lmax = 20;
reklas = 1;
l = 1;


while (l<lmax) && reklas
    X1pom = [];
    X2pom = [];
    X3pom = [];
    X4pom = [];
    reklas = 0;

    for i=1:N1
        d1 = 0.5*(X1(:,i)-M1)'*S1^(-1)*(X1(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X1(:,i)-M2)'*S2^(-1)*(X1(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        d3 = 0.5*(X1(:,i)-M3)'*S3^(-1)*(X1(:,i)-M3) + 0.5*log(det(S3)) - 0.5*log(P3);
        d4 = 0.5*(X1(:,i)-M4)'*S4^(-1)*(X1(:,i)-M4) + 0.5*log(det(S4)) - 0.5*log(P4);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X1(:,i)];
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X1(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X1(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X1(:,i)];
            reklas = 1;
        end
    end
    for i=1:N2
        d1 = 0.5*(X2(:,i)-M1)'*S1^(-1)*(X2(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X2(:,i)-M2)'*S2^(-1)*(X2(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        d3 = 0.5*(X2(:,i)-M3)'*S3^(-1)*(X2(:,i)-M3) + 0.5*log(det(S3)) - 0.5*log(P3);
        d4 = 0.5*(X2(:,i)-M4)'*S4^(-1)*(X2(:,i)-M4) + 0.5*log(det(S4)) - 0.5*log(P4);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X2(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X2(:,i)];
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X2(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X2(:,i)];
            reklas = 1;
        end
    end
    for i=1:N3
        d1 = 0.5*(X3(:,i)-M1)'*S1^(-1)*(X3(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X3(:,i)-M2)'*S2^(-1)*(X3(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        d3 = 0.5*(X3(:,i)-M3)'*S3^(-1)*(X3(:,i)-M3) + 0.5*log(det(S3)) - 0.5*log(P3);
        d4 = 0.5*(X3(:,i)-M4)'*S4^(-1)*(X3(:,i)-M4) + 0.5*log(det(S4)) - 0.5*log(P4);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X3(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X3(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X3(:,i)];
        else
            X4pom = [X4pom X3(:,i)];
            reklas = 1;
        end
    end
    for i=1:N4
        d1 = 0.5*(X4(:,i)-M1)'*S1^(-1)*(X4(:,i)-M1) + 0.5*log(det(S1)) - 0.5*log(P1);
        d2 = 0.5*(X4(:,i)-M2)'*S2^(-1)*(X4(:,i)-M2) + 0.5*log(det(S2)) - 0.5*log(P2);
        d3 = 0.5*(X4(:,i)-M3)'*S3^(-1)*(X4(:,i)-M3) + 0.5*log(det(S3)) - 0.5*log(P3);
        d4 = 0.5*(X4(:,i)-M4)'*S4^(-1)*(X4(:,i)-M4) + 0.5*log(det(S4)) - 0.5*log(P4);
        if d1==min([d1 d2 d3 d4])
            X1pom = [X1pom X4(:,i)];
            reklas = 1;
        elseif d2==min([d1 d2 d3 d4])
            X2pom = [X2pom X4(:,i)];
            reklas = 1;
        elseif d3==min([d1 d2 d3 d4])
            X3pom = [X3pom X4(:,i)];
            reklas = 1;
        else
            X4pom = [X4pom X4(:,i)];
        end
    end
    clear X1 X2 X3 X4
    X1 = X1pom;
    X2 = X2pom;
    X3 = X3pom;
    X4 = X4pom;
    
    %{

    figure();
    plot(X1(1,:), X1(2,:), 'ro', X2(1,:), X2(2,:), 'bo', X3(1,:), X3(2,:), 'go', X4(1,:), X4(2,:), 'yo');
    legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4');
    xlabel('x1');ylabel('x2');
    grid on;
    title(['Klasterizacija broj' num2str(l)]);
    pause
  
    %}

    N1 = length(X1);
    N2 = length(X2);
    N3 = length(X3);
    N4 = length(X4);
    P1 = N1/(2*N);
    P2 = N2/(2*N);
    P3 = N3/(2*N);
    P4 = N4/(2*N);
    M1 = mean(X1,2);
    M2 = mean(X2,2);
    M3 = mean(X3,2);
    M4 = mean(X4,2);
    S1 = cov(X1');
    S2 = cov(X2');
    S3 = cov(X3');
    S4 = cov(X4');
    l = l + 1;
end
ll(k) = l;
end
l_avg = mean(ll);
