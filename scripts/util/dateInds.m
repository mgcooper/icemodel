function [ starti, endi ] = dateInds(startdate,enddate,t,varargin)
%dateInds returns the start and end indices for the startdate and enddate
%on the calendar t. Note legacy behavior accepted startdate and enddate
%format yri,moi,dayi,hri, new behavior requires these be 1xn vectors e.g.
%[yri moi dayi hri]

% I tried to somplify this but i realized the reason it's complicated is
% becuase of variable input. by requiring the three inputs, we just have to
% put startdate and enddate in [] and it's much easier

% startdate can be a 1xn datenumber e.g. [yri moi dayi hri] or a 1x1
% datenum vector, same for enddate
% t can be a nx7 time object created by timebuilder or a nx1 datenum vector

% t = time_builder(1998,10,1,1999,9,30,24);

% if startdate and enddate are 1xn datenum vectors they are converted to
% datevec format

% tt = timebuilder(2015,7,1,1,0,2015,7,31,1,0,1);

%%%%%%%%% check inputs

if nargin == 9 % assume legacy format
    yri     =   startdate;
    moi     =   enddate;
    dayi    =   t;
    hri     =   varargin{1};
    yrf     =   varargin{2};
    mof     =   varargin{3};
    dayf    =   varargin{4};
    hrf     =   varargin{5};
    t       =   varargin{6};
    
    if size(t,2)==1 % assume a single vector of datenums
        t   =   datevec(t);
    end

    starti  =   find(t(:,1) == yri & t(:,2) == moi & t(:,3) == dayi &...
                    t(:,4) == hri, 1, 'first');
    endi    =   find(t(:,1) == yrf & t(:,2) == mof & t(:,3) == dayf &...
                    t(:,4) == hrf, 1, 'last');
else
    
    % must be 3 inputs
    narginchk(3,3)

    % convert startdate and enddate to rows
    if iscolumn(startdate) && numel(startdate)>1
        startdate = startdate';
    elseif iscolumn(enddate) && numel(enddate)>1
        enddate = enddate';
    end

    % startdate and enddate must be the same size
    assert(isequal(size(startdate),size(enddate)), ...
                'startdate must be the same size as enddate');

    % startdate/enddate can be vectors of datenums/datetimes or they can be
    % 1x7 datevectors for compatibility with timebuilder
    assert(size(startdate,2)==1 || (size(startdate,2)<7 && size(startdate,2)>2) , ...
            ['startdate and enddate must be 1x1 datenums, '             ...
            'or vectors of size 1xn where 2<n<7 of the form '           ...
            '[year, month, day, hour, minute, second] ']);
    % if start/enddate are datenums, convert to datetimes
    if size(startdate,2)==1 
        if ~isdatetime(startdate)
            assert(isdatetime(datetime(startdate,'ConvertFrom','datenum')), ...
                ['startdate and enddate must be vectors of size 1xn where 2<n<7 ' ...
                'of the form [year, month, day, hour, minute, second] ' ...
                'or 1x1 datenums, or 1x1 datetimes']);
            % if that worked, then this should work
            startdate   = roundn(startdate,-8); % round to nearsest ms
            startdate   = datetime(startdate,'ConvertFrom','datenum');
        end
    end
    % repeat for enddate
    if size(enddate,2)==1 
        if ~isdatetime(enddate)
            assert(isdatetime(datetime(enddate,'ConvertFrom','datenum')), ...
                ['startdate and enddate must be vectors of size 1xn where 2<n<7 ' ...
                'of the form [year, month, day, hour, minute, second] ' ...
                'or 1x1 datenums, or 1x1 datetimes']);
            % if that worked, then this should work
            enddate     = roundn(enddate,-8); % round to nearsest ms
            enddate     = datetime(enddate,'ConvertFrom','datenum');
        end
    end
    % repeat for t
    if size(t,2)==1 
        if ~isdatetime(t)
            assert(isdatetime(datetime(t,'ConvertFrom','datenum')), ...
                ['startdate and enddate must be vectors of size 1xn where 2<n<7 ' ...
                'of the form [year, month, day, hour, minute, second] ' ...
                'or 1x1 datenums, or 1x1 datetimes']);
            % if that worked, then this should work
            t   = roundn(t,-8); % round to nearsest ms
            t   = datetime(t,'ConvertFrom','datenum');
        end
    end
    % if nc>1, check if startdate, enddate, and t are datevecs and convert
    % to datetime
    nc                  =   size(t,2);
    switch nc
        case 7  % assume a timebuilder object, position 7 are datenums
            t           =   t(:,7);
    end % else assume a vector of datenums

    nc                  =   size(startdate,2);
    switch nc
        case {3, 6} % datevec format can have 3 or 7 columns
            startdate   =   datenum(startdate);
            enddate     =   datenum(enddate);
        case 4      % if [Y,M,D,H] was given, append the M,S
            startdate   =   datenum([startdate,0,0]);
            enddate     =   datenum([enddate,0,0]);
        case 5      % if [Y,M,D,H,M] was given, append the S
            startdate   =   datenum([startdate,0]);
            enddate     =   datenum([enddate,0]);
    end

    % round to nearest milisecond. I commented this out after adding the
    % conversion to datetimes in the checks above
%     if ~isdatetime(startdate); startdate = roundn(startdate,-8);end
%     if ~isdatetime(enddate); enddate = roundn(enddate,-8); end
%     if ~isdatetime(t); t = roundn(t,-8); end

    %%%%%%%%% find start and end indices
    starti  =   find(startdate == t, 1, 'first');
    endi    =   find(enddate >= t, 1, 'last');
end

end

