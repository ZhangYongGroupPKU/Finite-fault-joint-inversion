function [Xsolve,ob,syn,astftele,synastf] = shift_ob(tstele,tssm,Xsolve,ob,syn,obt,obsm,hdata_use,...
    meangreen,G_jt,sGDT,srate,ifastf)
   
    for k = 1:1
        if tstele
            syntele = syn(1:numel(obt));
            syntele = reshape(syntele,size(obt));
            [obstele,dttele,~] = rup_obssyn(reshape(hdata_use(1:numel(obt)),...
                size(obt)),syntele,srate);
            
            ob(1:numel(obstele)) = obstele*meangreen;
            hdata_use(1:numel(obt)) = obstele(:);
            if ifastf==1
            load('astfneed.mat');
            obt=reshape(obstele,size(gtele,1),size(gtele,4));
            
            [astftele,synastf]=getastftele(obt*meangreen,gtele,source,grida,fault,lensub,band_tele,s1);
            astftele=astftele/srate;
            astftele=astftele(:)/mobstf;
          
            hdata_use(end-size(sGDT,1)+1:end-size(sGDT,1)+size(astftele,1))=astftele;
            else 
                synastf=[];
                astftele=[];
            end
            %mean(dttele)
        else
            synastf=[];
            astftele=[];
        end
        if tssm
            synsm = syn(numel(obt)+1:end);
            synsm = reshape(synsm,size(obsm,1),[]);
            [obssm,dtsm,~] = rup_obssyn_sm(reshape(hdata_use(numel(obt)+1:...
                numel(ob)),size(obsm,1),[]),synsm,srate);
            ob(numel(obt)+1:end) = obssm*meangreen;
            hdata_use(numel(obt)+1:numel(ob)) = obssm(:);
            %mean(dtsm)
        end
        %save('dts.mat','dttele','dtsm');
        if tstele || tssm
            tic
            [Xsolve,~] = cgls_xu_invt0(G_jt,sGDT,hdata_use,1e-170,1000);
            toc
            syn = reshape(G_jt(1:numel(ob),:)*(Xsolve*meangreen),size(ob));
        else
            synastf=[];
            astftele=[];
        end
    end

end
