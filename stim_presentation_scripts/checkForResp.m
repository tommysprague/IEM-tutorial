function [resp, ts] = checkForResp(possResp, escapeKey)
% queries the keyboard to see if a legit response was made
% returns the response and the timestamp
resp = 0;
ts = nan;
[keyIsDown, secs, keyCode] = KbCheck;
if sum(keyCode)==1   % if at least one key was pressed
    keysPressed = find(keyCode);
    % in the case of multiple keypresses, just consider the first one
    if find(keysPressed(1)==possResp)
        resp = find(keyCode);
        ts = secs;
    end
    if keysPressed(1)==escapeKey
        Screen('CloseAll');
        resp = -1;
    end
end