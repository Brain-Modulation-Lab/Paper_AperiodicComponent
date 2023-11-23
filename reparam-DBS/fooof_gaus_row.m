function gaus_row = fooof_gaus_row(SUBJECT,SESSION_ID,electrode,type,fooofresult)
nr2=[];
gaus_row = table();
if ~isempty(fooofresult)
  ng = size(fooofresult.gaussian_params,1);
  nr2.subject=repmat({SUBJECT},ng,1);
  nr2.session_id=repmat(SESSION_ID,ng,1);
  nr2.electrode =repmat(electrode,ng,1);
  nr2.epoch_type=repmat(type,ng,1);
  nr2.gaus_id = (1:ng)';
  nr2.gaus_freq = fooofresult.gaussian_params(:,1);
  nr2.gaus_amp = fooofresult.gaussian_params(:,2);
  nr2.gaus_sigma = fooofresult.gaussian_params(:,3);  
  gaus_row = struct2table(nr2);
end
