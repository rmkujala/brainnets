function funpsy_interpvol(cfg);

    psess=cfg.psess;
    rois=psess.rois;
    blacklist=psess.blacklist;
    R=length(rois);
    vals=cfg.vals;

    vals(blacklist,:)=0;
    if(R ~= size(vals,1))
        error(strcat('mismatch in the length: vals: ', num2str(length(vals)), 'ROIS:', num2str(R)))
    end

	Nvals=size(vals,2);

	allvols=zeros([psess.datasize(1:3) Nvals]);




    for currval=1:Nvals

		vol=zeros(psess.datasize(1:3));
		%load(psess.groupmask)
		%vol(find(groupmask==1))=NaN;
		%compare=vol;
		X=zeros(R,1);
		Y=zeros(R,1);
		Z=zeros(R,1);
		Xi=[];
		Yi=[];
		Zi=[];
		for r=1:R
			c=rois(r).centroid;
			X(r)=c(1);
			Y(r)=c(2);
			Z(r)=c(3);


		    map=rois(r).map;
		    for m=1:size(map,1)
		        x=map(m,1);
		        y=map(m,2);
		        z=map(m,3);
                vol(x, y, z) = vals(r, currval);
            end
		%        compare(x,y,z)=vals(r);
		%        vol(x,y,z)=NaN;
		%        Xi=[Xi;x];
		%        Yi=[Yi;y];
		%        Zi=[Zi;z];
		%    end
			vol(c(1),c(2),c(3))=vals(r,currval);

        end

        [gX,gY,gZ] = meshgrid(2:3:109,2:3:91,2:3:91);
        [bgX, bgY,bgZ] = meshgrid(1:psess.datasize(2),1:psess.datasize(1),1:psess.datasize(3));
        voli=interp3(gX,gY,gZ,vol(2:3:91,2:3:109,2:3:91),bgX, bgY,bgZ,'nearest');

        voli(find(isnan(voli)))=0;

        amp=1;
%         if(cfg.box>0)
% gaussian leaves the original data a bit more intact...
%             voli=smooth3(voli,'gaussian',[ 3 3 3]);
%             %voli=smooth3(voli,'box',cfg.box);
%             % we need to scale the smoothed data
%             ids=find(vol>max([median(vals(:,currval)) 0]));
%             amp=mean(voli(ids)./vol(ids));
%             if(isnan(amp))
%                 save workspace_debug
%                 error('its nan!')
%             end
%         end

        allvols(:,:,:,currval)=voli/amp;

    end

    save_nii(make_nii(allvols),cfg.filename);


    %save_nii(make_nii(vol),cfg.filename);