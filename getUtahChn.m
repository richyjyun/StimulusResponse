%% Convert from Ripple channel to Utah channel
function chn = getUtahChn(c)

% Each connector on Ripple can have up to 128 channels
rem = mod((c-1),128)+1;
bank = (c-rem)/128;

% Get bank and channel within bank 
chn = bank*32 + rem;

end