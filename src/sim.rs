
//use math::round;

enum Field {
    UField,
    VField,
    SField
}

struct Fluid {
    density: f64,
    num_x: usize,
    num_y: usize,
    num_cells: i32,
    h: f64,            //size of cells
    u: Vec<f64>,      //velocity (horizontal)
    v: Vec<f64>,      //velocity (vertical)
    new_u: Vec<f64>,
    new_v: Vec<f64>,
    p: Vec<f64>,      //pressure
    s: Vec<i32>,      //scalar: 0: static object, 1: fluid
    m: Vec<f64>,      //smoke
    new_m: Vec<f64>
}

impl Fluid {
    fn new(density: f64, num_x: usize, num_y: usize, h: f64) -> Fluid {
        
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
                s: vec![0; num_cells], 
                m: vec![1.0; num_cells], 
                new_m: vec![0.0; num_cells] }

    }
    // Applys gravity to the fluid.
    // TODO:  Why skip when s is 0?
    fn integrate(&mut self, dt: f64, gravity: f64) {
        
        let n = self.num_y as usize;
        for i in 1..self.num_x as usize +1 {
            for j in 1..self.num_y as usize {
                if self.s[i*n + j] != 0 && self.s[i*n + j-1] != 0 { //if this cell, and the one below it are fluids, apply gravity
                    self.v[i*n + j] += gravity * dt;
                }
            }
        }
    }

    fn solve_incompressibility(&mut self, num_iters: i32, dt: f64) {
        let n = self.num_x;
        let cp = self.density * self.h / dt;

        for iter in 0..num_iters {
            for i in 1..self.num_x {
                for j in 1..self.num_y {
                    if self.s[i*n + j] == 0 {
                        continue; //if it is a static object...
                    }

                    //let s = self.s[i*n + j];
                    let sx0 = self.s[(i-1)*n + j];
                    let sx1 = self.s[(i+1)*n + j];
                    let sy0 = self.s[i*n + j-1];
                    let sy1 = self.s[i*n + j+1];
                    let s = sx0 + sx1 + sy0 + sy1;
                    
                    if s == 0 {
                        continue;
                    }

                    let div = self.u[(i+1)*n + j] - self.u[(i*n +j)] + self.v[i*n + j+1] - self.v[i*n + j];
                    let p = -div / s as f64;
                    //p *= scene.over_relaxation;
                    self.p[i*n + j] += cp * p;

                    self.u[i*n + j] -= sx0 as f64 * p;
                    self.u[(i+1)*n + j] += sx1 as f64 * p;
                    self.v[i*n + j] -= sy0 as f64 * p;
                    self.v[i*n + j+1] += sy1 as f64 * p;

                }

            }
        }
    }

    fn extrapolate(&mut self) {
        let n = self.num_y;
        for i in 0..self.num_x {
            self.u[i*n] = self.u[i*n + 1];
            self.u[i*n + self.num_y - 1] = self.u[i*n + self.num_y-2];
        }
        for j in 0..self.num_y {
            self.v[j] = self.v[n + j];
            self.v[(self.num_x+1)*n + j] = self.u[(self.num_x-2)*n +j]
        }
    }

    fn sample_field(&mut self, x: f64, y: usize, field: Field) {
        let n = self.num_y;
        let h = self.h;
        let h1 = 1.0 / h;
        let h2 = 0.5 * h;


        //let matcher = [x, self.num_x as f64 * &h ];
        let x = match [x, self.num_x as f64 * &h ].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        

        let x = match [x, h].iter().max_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let y = match [y, self.num_y * h as usize].iter().min() {Some(val) => val.clone(), None => 0};
        let y = match [y, h as usize].iter().max() {Some(val) => val, None => &0};

        let dx: f64 = 0.0;
        let dy: f64 = 0.0;

        //let x0 = 
        
        //match field {
        //    Field::UField => {let f = &mut self}
        //}


    }


}