clear all;
close all;
clc;

%% Generisanje odbiraka
N = 500;
%K1 je bimodalna Gauss-ova raspodela 0.6N1 + 0.4N2
%K2 je bimodalna Gauss-ova raspodela 0.6N3 + 0.4N4
M1 = [4 4]'; S1 = [2 -0.5;-0.5 2];
M2 = [4 8]'; S2 = [0.9 0.7;0.7 0.9];
M3 = [-5 6]'; S3 = [1.5 0.5;0.5 1.5];
M4 = [-5 2]'; S4 = [1 0.7;0.7 1];
X1 = zeros(N,2);
X2 = zeros(N,2);

for i=1:N
    x = rand();
    if x<0.6
        X1(i,:) = mvnrnd(M1,S1,1);
    else
        X1(i,:) = mvnrnd(M2,S2,1);
    end
end
for i=1:N
    x = rand();
    if x<0.6
        X2(i,:) = mvnrnd(M3,S3,1);
    else
        X2(i,:) = mvnrnd(M4,S4,1);
    end
end

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx');
legend('Klasa 1', 'Klasa 2');
xlabel('x1');ylabel('x2');
grid on;

%% Teorijska fgv

x_range = linspace(-10, 10, 100);
y_range = linspace(-4, 12, 100);
[X, Y] = meshgrid(x_range, y_range);

pdf1 = mvnpdf([X(:) Y(:)], M1', S1);
pdf2 = mvnpdf([X(:) Y(:)], M2', S2);
pdf_x1 = 0.6 * pdf1 + 0.4 * pdf2;


pdf3 = mvnpdf([X(:) Y(:)], M3', S3);
pdf4 = mvnpdf([X(:) Y(:)], M4', S4);
pdf_x2 = 0.6 * pdf3 + 0.4 * pdf4;

contour_levels_x1 = linspace(min(pdf_x1), max(pdf_x1), 10);
contour_levels_x2 = linspace(min(pdf_x2), max(pdf_x2), 10);

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx');
legend('Klasa 1', 'Klasa 2');
grid on;
hold all;
contour(X, Y, reshape(pdf_x1, size(X)), contour_levels_x1, 'r');
contour(X, Y, reshape(pdf_x2, size(X)), contour_levels_x2, 'b');
xlabel('x1'); ylabel('x2')


%% Bajesov klasifikator i procena gresaka

x = -10:0.1:10;
y = -4:0.1:12;
f1p = zeros(length(x),length(y)); 
f2p = f1p; 

c1 = 1/(2*pi*det(S1)^0.5);
c2 = 1/(2*pi*det(S2)^0.5);
c3 = 1/(2*pi*det(S3)^0.5);
c4 = 1/(2*pi*det(S4)^0.5);

e1 = 0;
e2 = 0;
for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i),y(j)]';
        f11 = c1*exp(-0.5*(X-M1)'*inv(S1)*(X-M1)); % Funkcije gustine verovatnoce
        f12 = c2*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
        f21 = c3*exp(-0.5*(X-M3)'*inv(S3)*(X-M3));
        f22 = c4*exp(-0.5*(X-M4)'*inv(S4)*(X-M4));
        f2p(i,j) = 0.6*f21 + 0.4*f22; 
        f1p(i,j) = 0.6*f11 + 0.4*f12;
        h = -log(f1p(i,j)/f2p(i,j));
        if h > 0
            e1 = e1 + f1p(i,j)*0.1^2;
        else
            e2 = e2 + f2p(i,j)*0.1^2;
        end
    end
end

h = -log(f1p./f2p); 

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx');
hold all;
contour(x,y,h',[0 0],'g','LineWidth',1.2);
legend('Klasa 1', 'Klasa 2');
grid on;

disp('Teorijska greška klasifikacije prve vreste: ');
disp(e1);
disp('Teorijska greška klasifikacije druge vreste: ');
disp(e2);
eu = e1+e2;
c = zeros(2,2);

for i = 1:N
    x = X1(i,:)';
    f11 = c1*exp(-0.5*(x-M1)'*inv(S1)*(x-M1));
    f12 = c2*exp(-0.5*(x-M2)'*inv(S2)*(x-M2));
    f21 = c3*exp(-0.5*(x-M3)'*inv(S3)*(x-M3));
    f22 = c4*exp(-0.5*(x-M4)'*inv(S4)*(x-M4));
    f2 = 0.6*f21 + 0.4*f22; 
    f1 = 0.6*f11 + 0.4*f12;
    h = -log(f1/f2);
    if h > 0
        c(1,2) = c(1,2) + 1;
    else
        c(1,1) = c(1,1) + 1;
    end
end
for i = 1:N
    x = X2(i,:)';
    f11 = c1*exp(-0.5*(x-M1)'*inv(S1)*(x-M1)); 
    f12 = c2*exp(-0.5*(x-M2)'*inv(S2)*(x-M2));
    f21 = c3*exp(-0.5*(x-M3)'*inv(S3)*(x-M3));
    f22 = c4*exp(-0.5*(x-M4)'*inv(S4)*(x-M4));
    f2 = 0.6*f21 + 0.4*f22; 
    f1 = 0.6*f11 + 0.4*f12;
    h = -log(f1/f2);
    if h > 0
        c(2,2) = c(2,2) + 1;
    else
        c(2,1) = c(2,1) + 1;
    end
end
disp(c);
%% Minimalna cena
c12 = 1;
c21 = 5;
x = -10:0.1:10;
y = -4:0.1:12;

h = -log(f1p./f2p); 
t = -log(c12/c21);
figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx');
hold all;
contour(x,y,h',[0 0],'m','LineWidth',1.2);
contour(x,y,h',[t t],'g','LineWidth',1.2);
legend('Klasa 1', 'Klasa 2','Bayes', 'Test minimalne cene');
grid on;

%% Neuman-Pearson
X = [X1; X2]';
Y_label = [ones(1,N), 2*ones(1,N)];
Y_pred = ones(1,2*N);

c1 = 1/(2*pi*det(S1)^0.5);
c2 = 1/(2*pi*det(S2)^0.5);
c3 = 1/(2*pi*det(S3)^0.5);
c4 = 1/(2*pi*det(S4)^0.5);

x = -10:0.1:10;
y = -4:0.1:12;
h = -log(f1p./f2p); 

cnt = 0;
for mi = 0.01:0.01:10
    cnt = cnt+1;
    eps2(cnt) = 0;
    for i = 1:length(x)
        for j = 1:length(y)
            if(h(i,j) <= -log(mi))
                eps2(cnt) = eps2(cnt) + 0.1*0.1*f2p(i,j);
            end
        end
    end
end

figure();
plot(0.01:0.01:10,eps2);
xlabel('\mu'); ylabel('\epsilon_2');
title('Zavisnost \epsilon_2 od \mu');


eps2_limit = 0.00003;
index = find(eps2 <= eps2_limit,1);
mi = 0.01:0.01:10;
mi = mi(index);

for i = 1:2*N
    xt = X(:,i);
    f11 = c1*exp(-0.5*(xt-M1)'*inv(S1)*(xt-M1)); 
    f12 = c2*exp(-0.5*(xt-M2)'*inv(S2)*(xt-M2));
    f21 = c3*exp(-0.5*(xt-M3)'*inv(S3)*(xt-M3));
    f22 = c4*exp(-0.5*(xt-M4)'*inv(S4)*(xt-M4));
    f2 = 0.6*f21 + 0.4*f22; 
    f1 = 0.6*f11 + 0.4*f12;
    h1 = -log(f1/f2);
    if h1 > -log(mi)
        Y_pred(i) = 2;
    end
end
    
figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx');
hold all;
contour(x,y,h',[0 0],'m','LineWidth',1.2);
contour(x,y,h',[-log(mi) -log(mi)],'g','LineWidth',1.2);
legend('Klasa 1', 'Klasa 2', 'Bayes', 'Neyman-Pearson');
grid on;

C = confusionmat(Y_label,Y_pred);
figure(10)
cm = confusionchart(Y_label,Y_pred);
e1 = C(1,2)/sum(C(1,:));
e2 = C(2,1)/sum(C(2,:));
disp('Greska prve vrste: ');
disp(e1);
disp('Greska druge vrste: ');
disp(e2);
%% Wald-ov sekvencijalni
e2 = 0.000000000001;
e1 = e2;

a = log(e2/(1-e1));
b = log((1 - e2)/e1);
Xp = [X1;X2];

N = 50;
figure()
title('Waldov sekvencijalni test');
hold all;
for m = 1:N
    
    Sm = 0;
    n1 = 0;
    Sm1 = [];
    while (Sm > a && Sm < b)
        
        ind = randi(1000);
        X = Xp(ind, :)';
        f11 = c1*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
        f12 = c2*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));% Funkcije gustine verovatnoce
        f21 = c3*exp(-0.5*(X-M3)'*inv(S3)*(X-M3));
        f22 = c4*exp(-0.5*(X-M4)'*inv(S4)*(X-M4));
        f1 = 0.6*f11 + 0.4*f12;
        f2 = 0.6*f21 + 0.4*f22;
        
        h = log(f2) - log(f1);
        Sm = Sm +h;
        n1 = n1 + 1;
        Sm1(n1) = Sm;
        
    end
    
    if Sm <= a
        plot(0:(n1), [0 Sm1], 'r')
        
    elseif Sm >= b
        plot(0:n1, [0 Sm1], 'b')
        
        e11 = e11 + 1;
    end
end
    
plot([0 10], [a a], 'm')
plot([0 10], [b b], 'c')


e1_niz = 1e-12:1e-8:1e-3;
greska1 = zeros(1,length(e1_niz));
for l = 1:length(e1_niz)
    e1 = e1_niz(l);
    e2 = e1;

    a = log(e2/(1-e1));
    b = log((1 - e2)/e1);

    N = 50;
    for m = 1:N    
        n1 = 0;
        Sm = 0;
        
        while (Sm > a && Sm < b)
            ind = randi(500);
            X = X1(ind, :)';
            f11 = c1*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
            f12 = c2*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));% Funkcije gustine verovatnoce
            f21 = c3*exp(-0.5*(X-M3)'*inv(S3)*(X-M3));
            f22 = c4*exp(-0.5*(X-M4)'*inv(S4)*(X-M4));
            f1 = 0.6*f11 + 0.4*f12;
            f2 = 0.6*f21 + 0.3*f22;

            h = log(f2) - log(f1);
            Sm = Sm +h;
            n1 = n1 + 1;
        end
        broj_odb1(m) = n1;
    end    
    greska1(l) = max(broj_odb1);
end
figure
semilogx(e1_niz, 10*log10(greska1) )
title('Zavisnost minimalnog broja odbiraka od greske tipa 1')

e2_niz = 1e-12:1e-8:1e-3;
greska2 = zeros(1,length(e2_niz));
for l = 1:length(e2_niz)
    e2 = e2_niz(l);
    e1 = e2;
    a = log(e2/(1-e1));
    b = log((1 - e2)/e1);
    N = 50;
    for m = 1:N 
        n2 = 0;
        Sm = 0;
        while (Sm > a && Sm < b)
            ind = randi(500);
            X = X2(ind, :)';
             f11 = c1*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
            f12 = c2*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));% Funkcije gustine verovatnoce
            f21 = c3*exp(-0.5*(X-M3)'*inv(S3)*(X-M3));
            f22 = c4*exp(-0.5*(X-M4)'*inv(S4)*(X-M4));
            f1 = 0.6*f11 + 0.4*f12;
            f2 = 0.6*f21 + 0.4*f22;

            h = log(f2) - log(f1);
            Sm = Sm +h;
            n2 = n2 + 1;
        end
        broj_odb2(m) = n2;
    end
    greska2(l) = max(broj_odb2);
end
figure
semilogx(e2_niz, 10*log10(greska2) )
title('Zavisnost minimalnog broja odbiraka od greske tipa 2')
