%% Coral Metapopulation Decline - Code Share
% Daniel M. Holstein - 2022
%
% This code runs an abundance-based metapopulation model for a coral species
% (Orbicella annularis) in the Western Atlantic, associated with the 
% manuscript "Predicting coral metapopulation decline in a changing 
% thermal environment". However, the code can easily be used generally.
%
% The goals of this model are to estimate persistence and fragmentation due 
% to thermal stress throughout the century, and to identify highly central
% nodes. Connecitivy information is derived from Holstein et al. (2014),
% thermal vulnerabilities are estimated from van Hooidonk et al. (2015).
% Habitat polygons and the percentage of reef in each polygon (Andrefouet &
% Guzman 2004), thermal vulnerabilities, and connectivity information
% are provided in LoadInputs.mat.

load LoadInputs.mat

% Not run
% fid=fopen('vulnerabilities');
% cdata=textscan(fid,'%f %f','delimiter',',');
% fclose(fid);
% 
% S = shaperead('8k_Polys.shp'); % may be in another dir
% Reefper = extractfield(S,'REEF_PER');
% 
% % Remove self-recruitment
% settle_noself = settle-diag(diag(settle));
% 
% % Remove nodes that experience no incoming or outgoing connections
% % ("baddies")
% baddies_settle = sum(settle_noself,1) == 0;
% baddies_source = sum(settle_noself,2) == 0;
% baddies = baddies_settle & baddies_source';
% baddies = find(baddies);
% ids = 1:3202;
% newids = ids;
% newids(baddies) = [];
% Reefper = Reefper(newids);
% vuln_good = vuln(newids,:);
%
% Row-normalize connectivity matrix
% settle_new = settle(newids,newids);
% Dmean = zeros(size(settle_new));
% for i = 1:length(newids)
%     for j = 1:length(newids)
%         Dmean(i,j) = settle_new(i,j)/sum(settle_new(i,:));
%     end
% end
% 
% Dmean(isnan(Dmean)) = 0;

% Calculate eigenvalue
DE = eigs(Dmean,1);

%% Run metapopulation model loop

% Calculate the alpha, or slope at the origin of the stock-recruitment curve
% at low population density (sensu White et al. 2010)
FLEP = 0.1; % FLEP can be changed
a = 1/(FLEP*DE);

% Calculate C at time 0.
C = zeros(length(newids),length(newids),9);
for i = 1:length(newids)
    for j = 1:length(newids)
        C(i,j,1) = Dmean(i,j)*Reefper(i)*a;
    end
end

% Project C over each decade
for d = 2:9
    for i = 1:length(newids)
        for j = 1:length(newids)
            C(i,j,d) = Dmean(i,j)*Reefper(i)*a*(1-(vuln_good(j,d)/10))*(1-vuln_good(i,d)/10);
        end
    end
end

% Calculate eigenvalues in each decade
E = zeros(1,9);
for d = 1:9
    [v e] = eigs(C(:,:,d),1);
    E(d) = e; 
end

% % Plot eigenvalue over time
% plot(E); hold on
% line([0 9],[1 1],'color','r');
% axis([1 9 0 3])

%% Look for SCCs, calculate centrality (pagerank)
for d = 1:9
    G(d).sp = sparse(C(:,:,d));
    G(d).G = digraph(C(:,:,d)-diag(diag(C(:,:,d)))); %rmove self loops for centrality measure
    G(d).PR = centrality(G(d).G,'pagerank','Importance',G(d).G.Edges.Weight);
    [So, Co] = graphconncomp(G(d).sp);
    G(d).S = So;
    G(d).C = Co;
    [v ed] = eigs(C(:,:,d),1);
    G(d).eig = ed;
    G(d).vec = v;
end

% See if SCCs are persistent
for d = 1:9
    for k = 1:G(d).S
        asnmnt = G(d).C;
        incomp = asnmnt == k;
        G(d).SSC(k).k = k;
        
        % If an SCC is composed by more than 2 nodes, check eigenvalue.
        % Otherwise, fill with NaN
        if sum(incomp) > 2
            G(d).SSC(k).eig = eigs(C(incomp,incomp,d),1);
        else
            G(d).SSC(k).eig = NaN;
        end
    end
end

%% Look at persistence for various values of FLEP

% Vector of FLEPs. This can be changed.
per = [.01 .05 .1 .15 .2];

for g = 1:length(per)
    a = 1/(per(g)*DE);

    % Calculate C at time zero
    Ce = zeros(length(newids),length(newids),9);
    for i = 1:length(newids)
        for j = 1:length(newids)
            Ce(i,j,1) = Dmean(i,j)*Reefper(i)*a;
        end
    end
    
    % Project C for each decade
    for d = 2:9
        for i = 1:length(newids)
            for j = 1:length(newids)
                Ce(i,j,d) = Dmean(i,j)*Reefper(i)*a*(1-(vuln_good(j,d)/10))*(1-vuln_good(i,d)/10);
            end
        end
    end
   
    for d = 1:9
        [v e] = eigs(Ce(:,:,d),1);
        Ee(d) = e;
    end
    p(g) = plot(Ee,'--','Color','k'); hold on
end
plot(E,'b-');
line([0 9],[1 1],'color','r');
axis([1 9 0 10])
axis square
xlabel('Decade')
ylabel('Eig')

