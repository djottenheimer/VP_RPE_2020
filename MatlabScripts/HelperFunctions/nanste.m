function  datste=nanste(dat,DM)

datste=nanstd(dat,[],DM)./sqrt(sum(~isnan(dat),DM));
%datste=nanstd(dat,[],DM)./sqrt(size(dat,DM));