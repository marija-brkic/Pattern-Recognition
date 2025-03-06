clear all;
close all;
clc;

N = 500;
M1 = [-3 -4]';
M2 = [3,4]';
M3 = [7,-3]';
S1 = [1.2 0.5; 0.5 1.2];
S2 = [0.6 0.2; 0.2 0.6];
S3 = [1 -0.7; -0.7 1];
 
X1 = mvnrnd(M1,S1,N);
X2 = mvnrnd(M2,S2,N);
X3 = mvnrnd(M3,S3,N);

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
legend('Klasa 1', 'Klasa 2', 'Klasa 3');
xlabel('x1');ylabel('x2');
grid on;

%% Linearni klasifikator-metod resupstitucije-klasifikacija jedne klase

X12 = [X1; X2];
[s_opt, v0_opt, Neps_opt, M1, M2, S1, S2] = metod_resupstitucije(X12, X3);

v_opt = (s_opt*S1+(1-s_opt)*S2)^(-1)*(M2-M1);
x11 = -8:0.1:12;
x21 = -(v0_opt+v_opt(1)*x11)/v_opt(2);

%%
figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);

Ylabel = [ones(1,N), 2*ones(1,N), 3*ones(1,N)];
Ypred = ones(1,3*N);

X = [X12; X3];

for i =1 : 3*N
   x = X(i,:)';
   if (v_opt' * x + v0_opt) > 0
       Ypred(i) = 3;  
   end  
end

ind = find(Ypred==1);
C = confusionmat(Ylabel,Ypred);

%% Klasifikacija druge dve klase

[s_opt, v0_opt, Neps_opt, M1, M2, S1, S2] = metod_resupstitucije(X1, X2);

v_opt = (s_opt*S1+(1-s_opt)*S2)^(-1)*(M2-M1);
x12 = -8:0.1:12;
x22 = -(v0_opt+v_opt(1)*x11)/v_opt(2);

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
plot(x12,x22,'k--');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija 1', 'Klasifikaciona linija 2');


if ~isempty(ind)
    X = X(ind,:);
    for i =1 :length(X)
        x = X(i,:)';
        if (v_opt' * x + v0_opt) > 0
            Ypred(i) = 2;  
        end  
    end
end

%% 
C = confusionmat(Ylabel,Ypred);
figure();
cm = confusionchart(Ylabel,Ypred);

%% Linearni klasifikator na bazi zeljenog izlaza

X12 = [X1; X2];

U = [-1*ones(1,2*N) ones(1,N); -X12', X3'];
G = ones(3*N,1); 
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3);

x11 = -8:0.1:12;
x21 = -(v0+v1*x11)/v2;

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);

Ylabel = [ones(1,N), 2*ones(1,N), 3*ones(1,N)];
Ypred = ones(1,3*N);

X = [X12; X3];
v_opt = [v1; v2];

v_opt = [v1; v2];
for i =1 : 3*N
   x = X(i,:)';
   if (v_opt' * x + v0) > 0
       Ypred(i) = 3;  
   end  
end

ind = find(Ypred==1);


U = [-1*ones(1,N) ones(1,N); -X1', X2'];
G = ones(2*N,1); 
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3);

x12 = -8:0.1:12;
x22 = -(v0+v1*x11)/v2;

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
plot(x12,x22,'k--');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija 1', 'Klasifikaciona linija 2');

v_opt = [v1; v2];
if ~isempty(ind)
    X = X(ind,:);
    for i =1 :length(X)
        x = X(i,:)';
        if (v_opt' * x + v0) > 0
            Ypred(i) = 2;  
        end  
    end
end

C = confusionmat(Ylabel,Ypred);
figure()
cm = confusionchart(Ylabel,Ypred);

%% Vece gamma

X12 = [X1; X2];

U = [-1*ones(1,2*N) ones(1,N); -X12', X3'];
G = [ones(2*N,1); 1.5 * ones(N,1)]; 
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3);

x11 = -8:0.1:12;
x21 = -(v0+v1*x11)/v2;

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);

Ylabel = [ones(1,N), 2*ones(1,N), 3*ones(1,N)];
Ypred = ones(1,3*N);

X = [X12; X3];
v_opt = [v1; v2];

v_opt = [v1; v2];
for i =1 : 3*N
   x = X(i,:)';
   if (v_opt' * x + v0) > 0
       Ypred(i) = 3;  
   end  
end

ind = find(Ypred==1);


U = [-1*ones(1,N) ones(1,N); -X1', X2'];
G = ones(2*N,1); 
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3);

x12 = -8:0.1:12;
x22 = -(v0+v1*x11)/v2;

figure();
plot(X1(:,1), X1(:,2), 'ro', X2(:,1), X2(:,2), 'bx', X3(:,1), X3(:,2), 'gh');
hold all;
plot(x11,x21,'k--');
plot(x12,x22,'k--');
xlabel('x1');ylabel('x2');
title('Klasifikacija');
grid on;
xlim([-8 12]);
ylim([-8 8]);
legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasifikaciona linija 1', 'Klasifikaciona linija 2');

v_opt = [v1; v2];
if ~isempty(ind)
    X = X(ind,:);
    for i =1 :length(X)
        x = X(i,:)';
        if (v_opt' * x + v0) > 0
            Ypred(i) = 2;  
        end  
    end
end

C = confusionmat(Ylabel,Ypred);
figure(10)
cm = confusionchart(Ylabel,Ypred);

%% Kvadratni klasifikator

N = 500;
X = zeros(2,N);
Y = zeros(2,N);
ro_1   = rand(1,N);
teta_1 = rand(1,N)*2*pi;
ro_2   = rand(1,N)+3;
teta_2 = rand(1,N)*2*pi;

X(1,:) = ro_1.*cos(teta_1);
X(2,:) = ro_1.*sin(teta_1);

Y(1,:) = ro_2.*cos(teta_2);
Y(2,:) = ro_2.*sin(teta_2);

figure();
hold all;
plot(X(1,:), X(2,:), 'ro', Y(1,:), Y(2,:), 'bx');
xlabel('x1');ylabel('x2');
grid on;
legend('Klasa 1', 'Klasa 2');

U = [-1*ones(1,N) ones(1,N); -X, Y; -X(1,:).^2 Y(1,:).^2; ...
    -X(2,:).^2 Y(2,:).^2; -2*X(1,:).*X(2,:) 2*Y(1,:).*Y(2,:)];
G = ones(2*N,1); %gamma
W = (U*U')\U*G;

v0 = W(1); v1 = W(2); v2 = W(3); v3 = W(4); v4 = W(5); v5 = W(6);

x1 = -5:0.1:5;
x2 = -5:0.1:5;
h = zeros(length(x1), length(x2));

for j = 1:length(x2)
    for i = 1:length(x1)
        h(i,j) = v0+v1*x1(i)+v2*x2(j)+v3*x1(i).^2+v4*x2(j).^2+2*v5*x2(j).*x1(i);
    end
end

figure();
hold all;
plot(X(1,:), X(2,:), 'ro', Y(1,:), Y(2,:), 'bx');
xlabel('x1');ylabel('x2');
grid on;
contour(x1,x2,h,[0 0], 'k--');
legend('Klasa 1', 'Klasa 2', 'Klasifikaciona linija');
title('Klasifikacija')