
//use math::round;

#[derive(PartialEq)]
enum Field {
    UField,
    VField,
    SField
}

pub struct Scene {
    pub gravity: f64,
    pub dt: f64,
    pub num_iters: i32,
    pub frame_nr: i32,
    pub over_relaxation: f64,
    pub obstacle_x: f64,
    pub obstacle_y: f64,
    pub obstacle_radius: f64,
    pub paused: bool,
    pub scene_nr: i32,
    pub show_obstacle: bool,
    pub show_stream_lines: bool,
    pub show_velocities: bool,
    pub show_pressure: bool,
    pub show_smoke: bool,
    pub sim_height: f64,
    pub fluid: Fluid
}

impl Default for Scene {
    fn default() -> Self {
        let mut sc: Scene =    Scene { gravity: -9.81, 
                dt: 1.0/60.0, 
                num_iters: 100, 
                frame_nr: 1, 
                over_relaxation: 1.9, 
                obstacle_x: 100.0, 
                obstacle_y: 100.0, 
                obstacle_radius: 0.15, 
                paused: false, 
                scene_nr: 0, 
                show_obstacle: false, 
                show_stream_lines: false, 
                show_velocities: false, 
                show_pressure: false, 
                show_smoke: false ,
                sim_height: 1.1,
                fluid: Fluid::new(100.0, 600, 600, 70.0)
            };
        
        sc.setup();
        return sc;
            
    }


}

impl Scene {
    pub fn setup(&mut self) {
        self.obstacle_radius = 0.15;
        self.over_relaxation = 1.9;
        let res: i32 = 100;
        let domain_height: f64 = 1.0;
        //let domain_width: f64 = domain_height / self.sim_height
        self.show_smoke = true;

        let n = self.fluid.num_y;

        for i in 0..self.fluid.num_x {
            for j in 0..self.fluid.num_y {
                let mut s = 1.0;
                if i==0 || j==0 || j==self.fluid.num_y-1 {
                    s = 0.0;
                }
                self.fluid.s[i*n + j] = s;
                if i==1 {
                    self.fluid.u[i*n +j] = 2.0;
                }


            }


        }

        let pipe_h = 0.1 * self.fluid.num_y as f64;
        let minj: f64 = 0.5 * self.fluid.num_y as f64 - 0.5*pipe_h;
        let minj: i32 = minj.floor() as i32;
        let maxj: f64 = 0.5 * self.fluid.num_x as f64 + 0.5*pipe_h;
        let maxj: i32 = maxj.floor() as i32;

       

        

        self.set_obstacle(6000.0, 12000.0);
        //self.set_square();
        //let sum: f64 = self.fluid.s.iter().sum();
        //println!("{:?}", sum);
        //for j in minj..maxj {
        //    self.fluid.u[(j + 500) as usize] = 1000.0;
        //}

        for i in 0..self.fluid.num_x {
            self.fluid.u[i + self.fluid.num_x as usize] = 100000.0;
        }



    }

    fn set_obstacle(&mut self, x: f64, y: f64,) {
        println!("setup obs");
        let vx: f64 = (x - self.obstacle_x) / self.dt;

        let vy: f64 = (y - self.obstacle_y) / self.dt;

        self.obstacle_x = x;
        self.obstacle_y = y;

        //let r = self.obstacle_radius;
        let r: f64 = 1000.0;

        let n = self.fluid.num_y;
        let cd = 2.0_f64.sqrt() * self.fluid.h;

        for i in 1..self.fluid.num_x-2 {
            for j in 1..self.fluid.num_y-2 {
                self.fluid.s[i*n + j] = 1.0;

                let dx = (i as f64 + 0.5) * self.fluid.h - x;
                let dy = (j as f64 + 0.5) * self.fluid.h - y;

                //println!("{}|{}", dx * dx + dy * dy, r * r);

                if dx * dx + dy * dy < r * r {
                    //println!("{}", "block");
                    self.fluid.s[i*n +j] = 0.0;
                    self.fluid.m[i*n +j] = 1.0;
                    //self.fluid.u[i*n +j] = 100.0;
                    //self.fluid.u[i*n +j] = 100.0;
                    //self.fluid.v[i*n +j] = 100.0;
                    //self.fluid.v[i*n +j+1] = 100.0;
                    //self.fluid.u[(i+1)*n + j] = vx / 2.0; //This one
                    //println!("vx={}", vx);
                    //self.fluid.v[i*n + j] = vy / 50.0;
                    //self.fluid.v[i*n + j+1] = vy / 100.0;

                }

            }
        }

    }

    fn set_square(&mut self) {
        let n = self.fluid.num_y;
        for i in 0..500 {
            for j in 0..500 {
                if i < 300 && i > 200 && j < 300 && j > 200 {
                    self.fluid.s[i*n + j] = 0.0;
                }
            }
        }
    }
}


pub struct Fluid {
    density: f64,
    pub num_x: usize,
    pub num_y: usize,
    num_cells: i32,
    h: f64,            //size of cells
    u: Vec<f64>,      //velocity (horizontal)
    v: Vec<f64>,      //velocity (vertical)
    new_u: Vec<f64>,
    new_v: Vec<f64>,
    p: Vec<f64>,      //pressure
    s: Vec<f64>,      //scalar: 0: static object, 1: fluid
    m: Vec<f64>,      //smoke
    new_m: Vec<f64>,
    //check: Vec<f64>,
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
                s: vec![0.5; num_cells], 
                m: vec![1.0; num_cells], 
                new_m: vec![0.0; num_cells],
            }

    }
    // Applys gravity to the fluid.
    // TODO:  Why skip when s is 0?
    fn integrate(&mut self, dt: f64, gravity: f64) {
        
        let n = self.num_y as usize;
        for i in 1..self.num_x - 1 as usize +1 {
            for j in 1..self.num_y -1 as usize {
                if self.s[i*n + j] != 0.0 && self.s[i*n + j-1] != 0.0 { //if this cell, and the one below it are fluids, apply gravity
                    self.v[i*n + j] += gravity * dt;
                }
            }
        }
    }

    fn solve_incompressibility(&mut self, num_iters: i32, dt: f64) {
        let n = self.num_x;
        let cp = self.density * self.h / dt;

        for _iter in 0..num_iters {
            for i in 1..self.num_x -1 {
                for j in 1..self.num_y {
                    if self.s[i*n + j] == 0.0 {
                        continue; //if it is a static object...
                    }

                    //let s = self.s[i*n + j];
                    let sx0 = self.s[(i-1)*n + j];
                    let sx1 = self.s[(i+1)*n + j];
                    let sy0 = self.s[i*n + j-1];
                    let sy1 = self.s[i*n + j+1];
                    let s = sx0 + sx1 + sy0 + sy1;
                    
                    if s == 0.0 {
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
        for j in 0..self.num_y - 1 {
            self.v[j] = self.v[n + j];
            self.v[(self.num_x-1)*n + j] = self.u[(self.num_x-2)*n +j]
        }
    }

    fn sample_field(&mut self, x: f64, y: f64, field: Field) -> f64 {
        let n = self.num_y;
        let h = self.h;
        let h1 = 1.0 / h;
        let h2 = 0.5 * h;


        let x = match [x, self.num_x as f64 * h ].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let x: f64 = match [x, h].iter().max_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let y = match [y, self.num_y as f64 * h ].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let y = match [y, h].iter().max_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        
        let mut dx: f64 = 0.0;
        let mut dy: f64 = 0.0;

        let f = match field {
            Field::UField => &self.u,
            Field::VField => &self.v,
            Field::SField => &self.s};


        match field {
            Field::UField=> dy = h2,
            Field::VField=> dx = h2,
            Field::SField=> {dx=h2; dy=h2},
            
        }
        
        let x0 = match [((x - dx)*h1).floor(), (self.num_x - 1) as f64 ].iter().min_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let tx = ((x-dx) - x0*h) *h1;
        let x1 = match [x0 + 1.0, self.num_x as f64 - 1.0].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let y0 = match [((y - dy)*h1).floor(), (self.num_y - 1) as f64 ].iter().min_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let ty = ((y-dy) - y0*h) *h1;
        let y1 = match [y0 + 1.0, self.num_y as f64 - 1.0].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};



        let sx = 1.0 - tx;
        let sy: f64 = 1.0 - ty;

        let val: f64 = 
        sx*sy * f[(x0 as usize * n + y0 as usize)]
        + tx*sy * f[x1 as usize*n + y0 as usize]
        + tx*ty * f[x1 as usize*n + y1 as usize]
        + sx*ty * f[x0 as usize*n + y1 as usize];

        //println!("{}",val);
        return val;


        
        //match field {
        //    Field::UField => {let f = &mut self}
        //}


    }


    fn avg_u(&self, i: usize, j: usize) -> f64 {
        let n = self.num_y;
        let u = (self.u[(i-1)*n + j] + self.u[i*n + j] + self.u[(i-1)*n + j+1] + self.u[i*n + j+1]) * 0.25;
        return u;
    }

    fn avg_v(&self, i: usize, j:usize) -> f64 {
        let n = self.num_y;
        let v = (self.v[(i-1)*n + j] + self.v[i*n + j] + self.v[(i-1)*n + j+1] + self.v[i*n + j+1]) * 0.25;
        return v;
    }

    fn advect_vel(&mut self, dt: f64) {
        self.new_u = self.u.to_vec();
        self.new_v = self.v.to_vec();

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.75 * h;

        for i in 1..self.num_x {
            for j in 1..self.num_y {
                if self.s[i*n + j] != 0.0 && self.s[(i-1)*n + j] != 0.0 && j < self.num_y -1 {
                    let mut x = i as f64 * h;
                    let mut y = j as f64 * h + h2;
                    let mut u = self.u[i*n + j];
                    let v = self.avg_v(i, j);
                    x = x - dt*u;
                    y = y - dt*v;
                    u = self.sample_field(x, y, Field::UField);
                    //print!("({})", u);
                    self.new_u[i*n + j] = u;
 
                }

                if self.s[i*n + j] != 0.0 && self.s[i*n + j-1] != 0.0 && i < self.num_x - 1 {
                    let mut x = i as f64 * h;
                    let mut y = j as f64 * h + h2;
                    let u = self.avg_u(i, j);
                    let mut v = self.v[i*n + j];
                    x = x - dt*u;
                    y = y - dt*v;
                    v = self.sample_field(x, y, Field::VField);
                    //print!("({})", v);
                    self.new_v[i*n + j] = v;
                } 
            }
        }

        self.u = self.new_u.to_vec();
        self.v = self.new_v.to_vec();
        //println!("{}", self.s[200]);

    }

    fn advect_smoke(&mut self, dt: f64) {
        //println!("running smoke");
       self.new_m = self.m.to_vec();
       let n = self.num_y;       
       let h = self.h;
       let h2 = 0.5 * h;

       for i in 1..self.num_x - 1 {
        for j in 1..self.num_y - 1 {
            if self.s[i*n + j] != 0.0 {
                let u = (self.u[i*n +j] + self.u[(i+1)*n + j]) * 0.5;
                let v = (self.v[i*n +j] + self.v[i*n + j+1]) * 0.5;
                let x: f64 = i as f64 * h + h2 - dt*u;
                let y: f64 = j as f64 * h + h2 - dt*v;

                self.new_m[i*n + j] = self.sample_field(x, y, Field::SField);
                //println!("{}", self.new_m[i*n + j]);


            }
        }
       }
       self.m = self.new_m.to_vec();
       //println!("{}|{}", self.m[250], self.m[2500]);
       
    }

    pub fn simulate(&mut self, dt:f64, gravity: f64, num_iters: i32) {
        //self.integrate(dt, gravity);
        //self.p = vec![0.0; self.num_cells as usize];
        self.solve_incompressibility(num_iters, dt);
        self.extrapolate();
        self.advect_vel(dt);
        self.advect_smoke(dt);
    }

    pub fn build_image(&self) -> Vec<u8> {
        let mut ctr: i32 = 0;
        let mut img: Vec<u8> = vec![];
        for i in 0..self.num_cells{
            //ctr += 1.0/250000.0;
        //for s in self.u.iter() {
            let mut color: f64 = 1.0 * (self.u[i as usize] + self.v[i as usize]) / 500.0;
            //if color < 0.0 {
            //    println!("neg");
            //}
            //if color > 255.0 {
            //    println!("Overflow");
            //}
            //if color > 10.0 {
            //    println!("{}", color);
            //}
            //let mut color: f64 = 255.0 * self.m[i as usize];
            //println!("{}", self.m[i as usize]);
            //let mut color: f64 = 0.000000005 * self.p[i as usize] + 100.0;
            //println!("{}", color);
            //let mut color: f64 = 255.0 * self.v[i as usize];

            /*let mut color: f64 = 255.0 * ctr;
            if i % 500 == 0 {
                
                color = 0.0;
            } else {
                println!("{}", i % 500);
                //color = 128.0;
                color = (i % 500) as f64;
            }
             */
            //println!("{}", color);
            color = color.floor();
            if self.s[i as usize] == 1.0 {
                
            img.push(color as u8);
            img.push(color as u8);
            img.push(color as u8);
            } else {
                img.push(255);
                img.push(0);
                img.push(0);
            }
            //println!("{}", self.s[i as usize]);
            
            
            //img.push(255);
        }

        let sum: f64 = self.s.iter().sum();
        //println!("{}", ctr);
        return img;
    }

}