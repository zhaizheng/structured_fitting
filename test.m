
c = -1:0.1:1;
y = c.^2/10;
X = [cos(c);sin(c)]+0.1*rand(2,length(c));
subplot(1,2,1)
plot(X(1,:),X(2,:),'*');
[Q, x0, Theta, Phi, re] = Factorization(X);
hold on
phi = min(Phi):(max(Phi)-min(Phi))/20:max(Phi);
[psi, M] = Psi(phi, Theta);
Z = x0 + Q*M;
plot(Z(1,:),Z(2,:),'-*');
subplot(1,2,2)
plot(re)



function [Q, x0, Theta, Tau, re] = Factorization(X)
    d = 1;
    D = size(X,1);
    CX = X-mean(X,2);
    [U,~,~] = svd(CX);
    Tau = U(:,1:d)'*CX;
    n = 1;
    Theta = zeros(d*(d+1)/2, D-d);
    Q = eye(D);
    while true
        [Q, x0, Theta] = Regression(Tau, X, Theta, Q);
        Tau_new = Projection(x0, Q, Theta, X, d, Tau);
        if norm(Tau_new-Tau)<1.e-5 || n > 100
            break;
        end
        Tau = Tau_new;
        [~, M] = Psi(Tau, Theta);
        re(n) = norm(X-x0-Q*M);
        n = n+1;
    end
end



function Tau = Projection(x0, Q, Theta, X, d, Tau) 
    n = size(X,2);
    D = size(X,1);
    U = Q(:,1:d);
    Up = Q(:,d+1:D);
    A = toTensor(Theta, d);
    for i = 1:n
        s = U'*(X(:,i)-x0);
        c = Up'*(X(:,i)-x0);
        tau = Tau(:,i);
        step = 5;
        while true 
            % while f_value(tau-(step*gradient(tau, s, c, A)), s, c, A) > f_value(tau,s,c,A)
            %     step = step/2;
            % end
            %norm(gradient(tau, s, c, A))
            [g,H ] = gradient(tau, s, c, A);
            %tau = tau - step*gradient(tau, s, c, A);
            tau = tau - H\g;
            if norm(gradient(tau, s, c, A))<1.e-7
                break;
            end
            %fprintf('projection tau %d,step: %f,norm of gradient%f\n', i, step, norm(gradient(tau, s, c, A)));
        end
        Tau(:,i) = tau;
        %fprintf('projection tau %d,norm of gradient%f\n', i, norm(gradient(tau, s, c, A)));
    end
end


function [g,H] = gradient(tau, s, c, A)
    dim = size(A,1);
    d = length(tau);
    r = zeros(dim,1);
    for i = 1:dim
        r(i) = tau'*squeeze(A(i,:,:))*tau;
    end
    Mr = zeros(dim,d);
    for i = 1:dim
        Mr(i,:)  = squeeze(A(i,:,:))*tau;
    end
    M2 = zeros(d,d);
    for i = 1:dim
        M2 = M2 + squeeze(A(i,:,:))*(r(i)-c(i));
    end
    g = 2*(tau-s)+4*Mr'*(r-c);
    H = 2*eye(d)+8*Mr'*Mr+4*M2;
end


function f = f_value(tau, s, c, A)
    dim = size(A,1);
    r = zeros(dim,1);
    for i = 1:dim
        r(i) = tau'*squeeze(A(i,:,:))*tau;
    end
    f = (s-tau)'*(s-tau)+(c-r)'*(c-r);
end


function A = toTensor(Theta, d)
    dh = size(Theta,2);
    A = zeros(dh,d,d);
    for i = 1:dh
        A(i,:,:) = vector2mat(Theta(:,i),d);
    end
end


function S = vector2mat(a, dim)
    M = zeros(dim);
    W = ones(dim);
    ind = triu(W)>0;
    M(ind) = a/2;
    S = M+M';
end


function [Q, c, Theta] = Regression(Phi, X, Theta, Q)
    d = size(Phi,1);
    D = size(X,1);
    while true
        [psi, M] = Psi(Phi, Theta);
        c = mean(X-Q*M,2);
        [U1, ~, U2] = svd((X-c)*M');
        Q_new = U1*U2';
        V = Q_new(:,d+1:D);
        Theta_new = (psi*psi')\psi*(X-c)'*V;
        error = norm(Theta_new-Theta)+norm(Q_new-Q);
        %fprintf('f_value%f,error=%f\n',norm(X-c-Q*M),error)
        if error < 1.e-8
            break;
        end
        Q = Q_new;
        Theta = Theta_new;
    end
end


function [psi, M] = Psi(Phi, Theta)
    d = size(Phi,1);
    n = size(Phi,2);
    psi = zeros(d*(d+1)/2,n);
    for i = 1:d
        for j = i:d
            psi((i-1)*d+j,:) = Phi(i,:).*Phi(j,:);
        end
    end
    M = [Phi; Theta'*psi];
end