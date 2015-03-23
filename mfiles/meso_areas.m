function [results, results_uf] = meso_areas(com1, com2, blacklistFileName, roisFileName)

load(roisFileName)
psess.rois=rois;
load(blacklistFileName)
psess.blacklist=blacklist;
%psess.groupmask='/Users/enrico/code/dyn_nets/stuff/mask.mat';
whitelist=setdiff(1:length(rois),blacklist);
coms = [com1; com2];
size(coms(1,:))


outstr={
    'best_NT'
    'best_AS'
    };
R=length(rois);


ids_to_atlasids=zeros(R,1);
atlas_size=zeros(10028,1);
for r=1:R
    atlas_ref{rois(r).atlas_ID}=rois(r).label;
    ids_to_atlasids(r)=[rois(r).atlas_ID];
    if(ismember(r,whitelist))
        % atlas_IDs range to 10028 apparently
        atlas_size(rois(r).atlas_ID)=   atlas_size(rois(r).atlas_ID)+1;
    end
end
        

dataALL=zeros(R,2);
for f=1:size(coms,1)
    com = coms(f,:);
    datatemp=com; %?
    data=zeros(R,1);
    data(whitelist)=double(datatemp); % added the plus 1 since 0 was a cluster
    dataALL(whitelist,f)=datatemp;
    
    ud=unique(datatemp);    %cluster labels
    C=length(ud); %number of communities
    if(f==1)
        udNT=ud;
        CNT=C;
    else
        udAS=ud;
        CAS=C;
    end
    for c=1:C
        ids=find(data==ud(c));
        b=1:length(atlas_ref);
        b=b';
        [h ]=histc(ids_to_atlasids(ids),b);
        disp(num2str(c))
        A=[h b];
        AA=flipud(sortrows(A,1));
        temp=AA(1,1);
        tempid=1;
        while temp>2
            disp([atlas_ref{AA(tempid,2)} ' ' num2str(AA(tempid,1)) ' ' num2str(AA(tempid,1)/atlas_size(AA(tempid,2)))])
            tempid=tempid+1;
            temp=AA(tempid,1);
        end
        %disp('%%%%%%')
    end
end

for r=1:CNT
    idNT=find(dataALL(:,1)==udNT(r));
    for c=1:CAS
        idAS=find(dataALL(:,2)==udAS(c));        
        ids=intersect(idNT,idAS);
        if(length(ids)<2) continue; end
        b=1:length(atlas_ref);
        b=b';
        h=histc(ids_to_atlasids(ids),b);
        disp(num2str(c))
        persize=h./atlas_size(b);
        persize(find(isnan(persize)))=0;
        A=[h b persize];
        AA=flipud(sortrows(A,3));
        temp=AA(1,1);
        tempid=1;
        tempstr=cell(1);
        tempstr2='';
        while temp>=1
            s  = s+AA(tempid,3);
            tempstr{tempid}=[atlas_ref{AA(tempid,2)} ' ' num2str(AA(tempid,1)) ' ' num2str(AA(tempid,3)) ' ' num2str(AA(tempid,2))];
            tempstr2=[tempstr2 atlas_ref{AA(tempid,2)} ', '];
            tempid=tempid+1;
            temp=AA(tempid,1);
        end

        results{r,c}=tempstr;
        results_uf{r,c}=tempstr2;
                clear tempstr
                clear tempstr2
                
        %disp('%%%%%%');
    end
end


