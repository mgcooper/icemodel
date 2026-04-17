function tt = rmttleapinds(tt)
   feb29 = month(tt.Time) == 2 & day(tt.Time) == 29;
   tt = tt(~feb29, :);
end
