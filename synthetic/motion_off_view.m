function y = motion_off_view(vf, u0, dt, y0)
% Simulate the motion of a Hill's vortex with a central velocity determined
% from the free-stream velocity 'u0', in the positive y direction, by
% removing velocity vectors in front per one frame. 'y0' specifies the
% current boundary of the vortex and 'dt' is the time step for motion.
% Central velocity.
uc = 3/2 * u0;
% y after decrement.
y = y0 - uc*dt;

% Remove front section.
if (y < vf.ybounds(1))
    idx = 0;
else
    idx = vf.getIndex_y(y);
end
vf.U_e(idx+1:end, :, :, :) = 0;
vf.updateFields();
