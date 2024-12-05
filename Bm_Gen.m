function Bm = Bm_Gen( M, T, S, SKm )

Bm=zeros(M,T);
for i=1:M
    for k=1:S
        Bm(i,SKm((i-1)*S+k))=1;
    end
end
Bm = sparse(Bm);

% j = 1;
% locx = zeros(M*S/2,1);
% locy = zeros(M*S/2,1);
% value = zeros(M*S/2,1);
% for i=1:M/2
%     for k=1:S
%         locx(j) = i;
%         locy(j) = SKm((i-1)*S+k)-floor(T/2);
%         value(j) = 1;
%         j =j+1;
%     end
% end
% Bm = sparse(locx, locy, value);

end