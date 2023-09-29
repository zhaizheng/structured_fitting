figure('Position', [1, 1, 800, 350])
x = -3:0.1:3;
y = -3:0.1:3;
[X, Y] = meshgrid(x,y);
Z = exp(-(X.^2+2*Y.^2));%*((X/5-Y/3).^2);
t = tiledlayout(2,3,"TileSpacing","compact");

nexttile
q = 1;
Z = exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3 ])
hT = title('q=1')
set(hT, 'FontSize', 14)
nexttile
q = 1/3;
Z = exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3 ])
hT = title('q=1/3')
set(hT, 'FontSize', 14)
nexttile
q = 1/30;
Z = exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
Y1 = zeros(size(x));
Z1 = exp(-(x.^2+2*Y1.^2)).^q;
plot3(x,Y1,Z1,'r-','Linewidth',2)
hold on
Y2 = [-3:0.1:-1/sqrt(8*q),1/sqrt(8*q):0.1:3];
X2 = zeros(size(Y2));
Z2 = exp(-(X2.^2+2*Y2.^2)).^q;
plot3(X2,Y2,Z2,'r-','Linewidth',2)
axis([-3 3 -3  3])
hT = title('q=1/30')
set(hT, 'FontSize', 14)
nexttile
q = -1/16;
Z = -exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
axis([-3 3 -3  3])
hT = title('q=-1/16')
set(hT, 'FontSize', 14)
nexttile
q = -1/8;
Z = -exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
axis([-3 3 -3  3])
hT = title('q=-1/8');
set(hT, 'FontSize', 14)

nexttile
q = -1;
Z = -exp(-(X.^2+2*Y.^2)).^q;
mesh(X,Y,Z)
hold on
x = -sqrt(-1/2/q):0.1:sqrt(-1/2/q);
Y1 = zeros(size(x));
Z1 = -(exp(-(x.^2))).^(q);
plot3(x,Y1,Z1,'r-','Linewidth',2)
axis([-3 3 -3  3])
hT = title('q=-1');
set(hT, 'FontSize', 14)
exportgraphics(t,'demo_simple.eps','Resolution',400)