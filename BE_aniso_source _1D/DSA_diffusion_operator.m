%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffusion_operator = DSA_diffusion_operator(...
                               dt,...
							   x,dx,...
							   N_x,poly_order,...
                               mat_mass,...
                               mat_total_cross_section,...
                               mat_scattering_cross_section...
							   )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global reference_gradient

% This solver actually solves the P1 equation
% \partial_x g + \sigma_a \rho = \sigma_s Source
% 1/3\partial_x \rho + \sigma_t g = 0

% cross sections
mat_absorption_cross_section = mat_total_cross_section-mat_scattering_cross_section;

% upwind operators for rho and the first order moment g
% assemble matrix for Dx with central flux and the jump term
mat_D_plus = assemble_D_plus(N_x,poly_order,reference_gradient);
mat_D_minus = assemble_D_minus(N_x,poly_order,reference_gradient);
mat_D_central = 0.5*(mat_D_plus+mat_D_minus);
mat_D_jump = mat_D_plus-mat_D_minus;

% contribution of the jump terms, use exact integrated value here,
% one can also consider using numerical integral determined by the reference integration


% mat_g_operator_rho = -0.25*mat_D_jump;
% mat_rho_operator_g = -0.375*mat_D_jump;
mat_g_operator_rho = -0.5*mat_D_jump;
mat_rho_operator_g = -0.75*mat_D_jump;
%mat_rho_operator_g = -1.125*mat_D_jump;

% diffusion operator
%diffusion_operator = -mat_D_central*1/3*inv(mat_mass/dt+mat_total_cross_section+mat_rho_operator_g)*mat_D_central...
%			         +mat_absorption_cross_section+mat_mass/dt...
%				     +mat_g_operator_rho;
%{
diffusion_operator = -mat_D_central*1/3*inv(eps*mat_mass/dt+mat_total_cross_section/eps)*mat_D_central...
			         +mat_absorption_cross_section+eps*mat_mass/dt...
				     +mat_g_operator_rho;
%}


diffusion_operator =-mat_D_central*1/3*inv(mat_mass/dt+mat_total_cross_section)*mat_D_central...
 			         +mat_absorption_cross_section+mat_mass/dt...
 				     +mat_g_operator_rho;

%diffusion_operator = 1/3 * (mat_D_central - 1/2 *mat_D_jump)*(-inv(mat_total_cross_section)*mat_D_central+(1/2)*inv(mat_total_cross_section)*mat_D_jump)+mat_absorption_cross_section;


%diffusion_operator = diffusion_operator;



% consistent
%{
diffusion_operator = -dt^2*mat_D_central*1/3*inv(mat_mass+dt*(mat_total_cross_section+mat_rho_operator_g))*mat_D_central...
			         +dt*(mat_absorption_cross_section+mat_g_operator_rho)...
				     +mat_mass;
%}

% partially consistent
%{
diffusion_operator = -dt^2*mat_D_central*1/3*inv(mat_mass+dt*mat_total_cross_section)*mat_D_central...
			         +dt*(mat_absorption_cross_section...
				         +mat_g_operator_rho)...
                      +mat_mass;
diffusion_operator = mat_mass
%}

%{
diffusion_operator = (-mat_D_central*1/3*inv(mat_total_cross_section+mat_mass/dt)*mat_D_central...
			          +mat_absorption_cross_section+mat_mass/dt...
				      +mat_g_operator_rho);
%}
end


