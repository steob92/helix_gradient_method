#[macro_use]
use peroxide::fuga::*;

pub const N_AIR : f64 = 1.0003;


// Equation for a 1D parabola
pub fn parabola(
    x : f64, x_shift : f64, p : &[f64]
) -> f64 {
    let x_adj = x - x_shift - 45.;
    p[2] * (x_adj).powi(2) + p[1] * (x_adj) + p[0]
}

// Equation of a 2D parabola
pub fn parabola_2d(
    x : f64, y : f64, x_shift : f64, y_shift : f64,
    p : &[f64]
) -> f64 {
    let xi = x - x_shift - 45.;
    let yi = y - y_shift - 45.;

    let mut z = p[0] + p[1]*xi + p[2]*xi*xi + p[3]*yi  ; 
    z += p[4]*yi*yi + p[5]*xi*yi +p[6]*xi*xi*yi + p[7]*xi*yi*yi +p[8]*xi*xi*yi*yi;
    z
}


pub fn params_from_2dx(
        x : f64,x_shift : f64, p : &[f64], 
    ) -> (f64, f64, f64){
        let xi = x - x_shift - 45.;
        let a0 = p[0] + p[1]*xi + p[2]*xi*xi;
        let a1 = p[3] + p[5]*xi + p[6]*xi*xi;
        let a2 = p[4] + p[7]*xi + p[8]*xi*xi;
        (a0, a1, a2)
    }

pub fn params_from_2dy( 
        y : f64, y_shift : f64, p :&[f64],
    ) -> (f64, f64, f64) {
        let yi = y - y_shift - 45.;
        let a0 = p[0] + p[3]*yi + p[4]*yi*yi;
        let a1 = p[1] + p[5]*yi + p[7]*yi*yi;
        let a2 = p[2] + p[6]*yi + p[8]*yi*yi;
        (a0, a1, a2)
    }



// Using Snells law to get the angle of refraction
pub fn snells_law(n1 : f64, n2 : f64, theta1 : f64) -> Option<f64> {
    match n2 == 0. {
        true => None,
        false => {
            Some((n1 * theta1.sin() / n2).asin())
        }
    }
}

pub fn line( x : f64, m : f64, c : f64) -> f64{
    m*x +c
}

// Use the dot product to get the angle between two vectors
pub fn get_angle(a : &Vec<f64>, b : &Vec<f64>) -> f64 {
    let mut norm_a = 0.;
    let mut norm_b = 0.;
    
    let mut dot = 0.;

    for i in 0..a.len() {
        dot+= a[i] * b[i];
        norm_a += a[i].powi(2);
        norm_b += b[i].powi(2);
    }
        
    (dot / norm_a.sqrt() / norm_b.sqrt()).acos()
}


// Get the normal vector of a generic parabola
pub fn get_norm(x : f64, p0 : f64, p1 : f64, p2 :f64) -> Vec<f64>{
    vec![2.*p2 * x + p1, -1.]
}
    

// Distance between a line and a parabola
// define the line in terms of the initial parameters
pub fn dist( x : f64, theta : f64, y0 : f64,  p : &[f64], rad_offset_x : f64) -> f64{
    let y_line = line(x, theta.tan(), y0);
    let y_curve = parabola(y_line, rad_offset_x, p);
    (x - y_curve).abs()

}
