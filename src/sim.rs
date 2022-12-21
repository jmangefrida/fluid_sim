struct Fluid {
    density: f64,
    num_x: i32,
    num_y: i32,
    num_cells: i32,
    h: f64,
    u: Vec<f64>,
    v: Vec<f64>,
    new_u: Vec<f64>,
    new_v: Vec<f64>,
    p: Vec<f64>,
    s: Vec<f64>,
    m: Vec<f64>,
    new_m: Vec<f64>
}

impl Fluid {
    fn new(density: f64, num_x: i32, num_y: i32, h: f64) -> Fluid {
        
        let num_cells: usize = num_x as usize * num_y as usize;
        
        Fluid { density: density,
                num_x: num_x, 
                num_y: num_y, 
                num_cells: num_cells as i32, 
                h: h, 
                u: vec![0.0; num_cells], 
                v: vec![0.0; num_cells], 
                new_u: vec![0.0; num_cells], 
                new_v: vec![0.0; num_cells], 
                p: vec![0.0; num_cells], 
                s: vec![0.0; num_cells], 
                m: vec![1.0; num_cells], 
                new_m: vec![0.0; num_cells] }

    }
    // Applys gravity to the fluid.
    // TODO:  Why skip when s is 0?
    fn integrate(&mut self, dt: f64, gravity: f64) {
        
        let n = self.num_y as usize;
        for i in 1..self.num_x as usize +1 {
            for j in 1..self.num_y as usize {
                if self.s[i*n + j] != 0.0 && self.s[i*n + j-1] != 0.0 {
                    self.v[i*n + j] += gravity * dt;
                }
            }
        }
    }
}