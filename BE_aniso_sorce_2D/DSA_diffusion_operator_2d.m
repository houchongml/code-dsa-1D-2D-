%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffusion_operator = DSA_diffusion_operator_2d(...
                               dt, ...
                               x, y, dx, dy, ...
                               N_x, N_y, ...
                               mat_total_cross_section, ...
                               mat_scattering_cross_section)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global reference_gradient_x reference_gradient_y

% This solver actually solves the P1 equation
% \partial_x g_x + \partial_y g_y + \sigma_a \rho = \sigma_s Source
% 1/3 \partial_x \rho + \sigma_t g_x = 0
% 1/3 \partial_y \rho + \sigma_t g_y = 0

% Cross sections
mat_absorption_cross_section = mat_total_cross_section - mat_scattering_cross_section;

% Upwind operators for rho and the first order moment g
% Assemble matrix for Dx and Dy with central flux and the jump term
mat_D_plus_x = assemble_D_plus_2d_x(N_x, N_y, reference_gradient_x);
mat_D_minus_x = assemble_D_minus_2d_x(N_x, N_y, reference_gradient_x);
mat_D_central_x = 0.5 * (mat_D_plus_x + mat_D_minus_x);
mat_D_jump_x = mat_D_plus_x - mat_D_minus_x;

mat_D_plus_y = assemble_D_plus_2d_y(N_x, N_y, reference_gradient_y);
mat_D_minus_y = assemble_D_minus_2d_y(N_x, N_y, reference_gradient_y);
mat_D_central_y = 0.5 * (mat_D_plus_y + mat_D_minus_y);
mat_D_jump_y = mat_D_plus_y - mat_D_minus_y;

% Contribution of the jump terms, use exact integrated value here,
% one can also consider using numerical integral determined by the reference integration
phi11 = -0.5* inv(-(2/pi)*mat_D_jump_x - (1\pi)*mat_D_jump_y +1.5*mat_total_cross_section) *mat_D_central_x;
phi12 = -0.5* inv(-(2/pi)*mat_D_jump_y - (1\pi)*mat_D_jump_x +1.5*mat_total_cross_section) *mat_D_central_y;
% Diffusion operator
diffusion_operator = 1.5*mat_D_central_x*phi11 - (1\pi)*mat_D_jump_x + 1.5*mat_D_central_y*phi12 - (1\pi)*mat_D_jump_y+mat_absorption_cross_section;

diffusion_operator = diffusion_operator;

% Consistent
%{
diffusion_operator = -dt^2 * (mat_D_central_x + mat_D_central_y) * (1 / 3) ...
                     * inv(mat_mass + dt * (mat_total_cross_section + mat_rho_operator_g)) ...
                     * (mat_D_central_x + mat_D_central_y) ...
                     + dt * (mat_absorption_cross_section + mat_g_operator_rho) ...
                     + mat_mass;
%}

% Partially consistent
%{
diffusion_operator = -dt^2 * (mat_D_central_x + mat_D_central_y) * (1 / 3) ...
                     * inv(mat_mass + dt * mat_total_cross_section) ...
                     * (mat_D_central_x + mat_D_central_y) ...
                     + dt * (mat_absorption_cross_section + mat_g_operator_rho) ...
                     + mat_mass;
diffusion_operator = mat_mass;
%}

% Alternative approach
%{
diffusion_operator = - (mat_D_central_x + mat_D_central_y) * (1 / 3) ...
                     * inv(mat_total_cross_section + mat_mass / dt) ...
                     * (mat_D_central_x + mat_D_central_y) ...
                     + mat_absorption_cross_section + mat_mass / dt ...
                     + mat_g_operator_rho;
%}
end