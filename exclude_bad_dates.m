% excluding bad jGCaMP7f data (earliest 20191104)
warning('NOTE: USING EXCLUDE_BAD_DATES, EXCLUDING 20191104')
idx_7f = find(ismember({mutant.construct}, '10.921'));
gcamp7f = mutant(idx_7f);
dates = gcamp7f.date;
dates_int = cellfun(@(x) str2double(x), dates);
bad_dates = dates_int >= 20191104; % exclude after this date
gcamp7f.df_fpeak_med = gcamp7f.df_fpeak_med(:,~bad_dates);
gcamp7f.decay_half_med = gcamp7f.decay_half_med(:,~bad_dates);
gcamp7f.rise_half_med = gcamp7f.rise_half_med(:,~bad_dates);
gcamp7f.timetopeak_med = gcamp7f.timetopeak_med(:,~bad_dates);
gcamp7f.decay_half_med_comp = gcamp7f.decay_half_med_comp(:,~bad_dates);
gcamp7f.rise_half_med_comp = gcamp7f.rise_half_med_comp(:,~bad_dates);
gcamp7f.f0 = gcamp7f.f0(~bad_dates);
gcamp7f.plate = gcamp7f.plate(~bad_dates);
gcamp7f.well = gcamp7f.well(~bad_dates);
gcamp7f.date = gcamp7f.date(~bad_dates);
gcamp7f.mCherry = gcamp7f.mCherry(~bad_dates);
gcamp7f.nSegment = gcamp7f.nSegment(~bad_dates);
gcamp7f.fmean = gcamp7f.fmean(~bad_dates);
gcamp7f.fmean_med = gcamp7f.fmean_med(:,:,~bad_dates);
gcamp7f.fmax_med = gcamp7f.fmax_med(~bad_dates);
gcamp7f.df_f_med = gcamp7f.df_f_med(:,:,~bad_dates);
gcamp7f.df_fpeak_med_comp = gcamp7f.df_fpeak_med_comp(:,~bad_dates);
gcamp7f.df_fnoise = gcamp7f.df_fnoise(:,~bad_dates);

mutant(idx_7f) = gcamp7f;