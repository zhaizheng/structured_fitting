test_diff_weight

%%
data = X2;
sigma=0.3;
%nexttile

figure('Position',[1,1,800,350])
t = tiledlayout(2,3,'TileSpacing','Compact');



%plot(X2(1,1:1000),X2(2,1:1000),'.','MarkerSize',10);

%hT = title('Original Data','interpreter','tex');
%set(hT, 'FontSize', 14)
%Title = {'log(x)','p^{1/4}(x)','p^{1/2}(x)','p^{3/4}(x)','p(x)'};
%Title = {'\Gamma(-10000,x)','\Gamma(-1,x)','\Gamma(0,x)','\Gamma(0.5,x)','\Gamma(1,x)'};
Title = {'q = -10^3','q=-10^2','q=-10','q=-1','q=0','q=1'};

for j = 1:6
    nexttile
    q = 1;
    [X, Y, Z] = draw1(data,sigma,q);
    mesh(X, Y, Z)
    %plot(T(1,:),T(2,:),'.');
    hold on
    z = [];
    for  i = 1: size(Re{7-j},2)
        z = [z, expfun(Re{7-j}(:,i),data,sigma,q)];
    end
    plot3(Re{7-j}(1,sign{7-j}<0.045),Re{7-j}(2,sign{7-j}<0.045),z(sign{7-j}<0.045),'r.','MarkerSize',5);
    hT = title(Title{7-j},'interpreter','tex');
    axis([-2, 2, -2 2]);
    %legend({'truth', 'ridge'})
    %set(gca,'fontsize',14)
    set(hT, 'FontSize', 14)
end

exportgraphics(t,'Nonlinear_Demo_3D.eps','Resolution',400)

figure('Position',[1,1,800,175])
t1 = tiledlayout(1,6,'TileSpacing','Compact');
for j = 1:6
    nexttile
    plot(Re{7-j}(1,sign{7-j}<0.045),Re{7-j}(2,sign{7-j}<0.045),'r.','MarkerSize',5);
    hT = title(Title{7-j},'interpreter','tex');
    axis([-2, 2, -2 2]);
end
exportgraphics(t1,'Nonlinear_Demo_2D.eps','Resolution',400)




function [X, Y, Z] = draw1(data,sigma,q)
    x = -2:0.1:2;
    y = -2:0.1:2;
    [X, Y] = meshgrid(x,y);
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            Z(i,j) = expfun([X(i,j);Y(i,j)], data, sigma,q);
        end
    end
end



function y = expfun(x, data, sigma, q)
    y = 0;
    for i = 1:size(data,2)
        y = y+ (1/3)/(sigma^3)* exp(-norm(x - data(:,i)).^2/(sigma^2));
    end
    y = y^q;
end

