clear 
figure('Position',[1,1,800,350])
s = {'1','1/3','1/30','-1/16','-1/8','-1'};
t = tiledlayout(2,3,'TileSpacing','Compact');
nexttile
rep(1,s{1})
nexttile
rep(1/3,s{2})
nexttile
rep(1/30,s{3})
nexttile
rep(-1/16,s{4})
nexttile
rep(-1/8,s{5})
nexttile
rep(-1,s{6})
exportgraphics(t,'two_peaks.eps','Resolution',400)

function rep(q,s)
a = [-1.5,0];
b = [1.5, 0];
x = -4:0.1:4;
y = -4:0.1:4;

[X,Y] = meshgrid(x,y);
Tra = cell(size(X,1),size(X,2));
for i = 1:size(X,1)
    for j = 1:size(X,2)  
        Z(i,j) = ff([X(i,j),Y(i,j)],a,q);%(f([X(i,j),Y(i,j)],a)+f([X(i,j),Y(i,j)],b))^q;
    end
end
mesh(X,Y,Z)
hold on

init_points = rand(2,6000);
init_points(1,:) = (init_points(1,:)-0.5)*8;
init_points(2,:) = (init_points(2,:)-0.5)*8;
ind = [];
for i = 1:6000
    [R(:,i), Tra{i}, e] = move(init_points(:,i), a, b, q);
    if e<0
        ind = [ind,i];
    end
    NEW_Z(i) = ff([R(1,i),R(2,i)],a,q);%(f([R(1,i),R(2,i)],a)+f([R(1,i),R(2,i)],b))^q;
end
hold on
%ind = R(1,:)~=0;
plot3(R(1,ind),R(2,ind),NEW_Z(ind),'r.')
hT = title(['q=',s]);
set(hT, 'FontSize', 14)
axis([-4 4 -4 4 ])
end
% hold on
% for i = 1:size(X,1)
%     Tra_value = zeros(1,size(Tra{i},2));
%     for j = 1:size(Tra{i},2)
%         Tra_value(j) = f([Tra{i}(1,j),Tra{i}(2,j)],a)+f([Tra{i}(1,j),Tra{i}(2,j)],b);
%     end
%     hold on
%     if length(Tra{i})>1
%         plot3(Tra{i}(1,:),Tra{i}(2,:),Tra_value(:),'.')
%     end
% end






function [x, X, e] = move(x,a,b,q)
     max_iter = 10^5;
     k = 1;
     X = [];
     while true
        g = gf2(x,a,b,q);
        [U,E] = eig(Hf2(x,a,b,q));
        [~, ind] = sort(diag(E),'descend');
        e = E(ind(2),ind(2));
        P = eye(length(x)) - U(:,ind(1))*U(:,ind(1))';
        x = x + 0.3*P * g;
        if norm(P*g)< 10^(-7)
            break;
            k
        end  
        k = k+1;
        X = [X,x];
     end
end



function s = Hf2(x,a,b,q)
    s = (Hf(x,a)+Hf(x,b))+ ((q-1)/(f(x,a)+f(x,b)))*(gf(x,a)+gf(x,b))*(gf(x,a)+gf(x,b))';
end

function g = gf2(x,a,b,q)
    g = (f(x,a)+f(x,b))^(-1)*(gf(x,a)+gf(x,b));
    % (f(x,a)+f(x,b))^(q-1)*(gf(x,a)+gf(x,b));
    %g = (gf(x,a)+gf(x,b));
end


function s = Hf(x,a)
    s = f(x,a)*([-2,0;0,-4] + ...
        [4*(x(1)+a(1))^2, 8*(x(1)+a(1))*(x(2)+a(2));  8*(x(1)+a(1))*(x(2)+a(2)), 16*(x(2)+a(2))^2]);
end

function s = gf(x,a)
    s = f(x,a)*[-2*(x(1)+a(1)),-4*(x(2)+a(2))]';
end

function s = ff(x,a,q)
    
    s = sign(q)*(f(x,a)+f(x,-a))^q;
end



function s = f(x,a)
    s = exp(-(x(1)+a(1))^2-2*(x(2)+a(2))^2);
end



