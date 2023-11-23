function aper_row = fooof_aper_row(SUBJECT,SESSION_ID,electrode,type,fooofresult)

nr=[];
nr.subject=SUBJECT;
nr.session_id=SESSION_ID;
nr.electrode = electrode;
nr.epoch_type=type;
if ~isempty(fooofresult)
  nr.aper_offset = fooofresult.aperiodic_params(1);
  nr.aper_knee = fooofresult.aperiodic_params(2);
  nr.aper_exp = fooofresult.aperiodic_params(3);
  nr.error = fooofresult.error;
  nr.r_squared = fooofresult.r_squared;
end
aper_row = struct2table(nr);