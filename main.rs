use std::{collections::vec_deque, f64::consts::PI, vec};

fn main() {
    let target = Matrix {
        rows: 2,
        cols: 1,
        value: vec![
            vec![10.0],
            vec![12.0]
        ]
    };
    let jacobian1 = jacob(8.0, 5.0, 3.0, 90.0*PI/180.0, 0.0, 0.0);
    ik_jpsd(jacobian1, 0.1, xy(11.0, 5.0, 2.0, 90.0*PI/180.0, 0.0, 0.0), target, 30);
}

fn jacob(l1: f64, l2: f64, l3: f64, t1: f64, t2: f64, t3: f64) -> Matrix {
    return Matrix {
        rows: 2,
        cols: 3,
        value: vec![
            vec![-l1*f64::sin(t1)-l2*f64::sin(t1+t2)-l3*f64::sin(t1+t2+t3),
                -l2*f64::sin(t1+t2)-l3*f64::sin(t1+t2+t3), -l3*f64::sin(t1+t2+t3)],
            vec![l1*f64::cos(t1)+l2*f64::cos(t1+t2)+l3*f64::cos(t1+t2+t3),
                l2*f64::cos(t1+t2)+l3*f64::cos(t1+t2+t3), l3*f64::cos(t1+t2+t3)]
        ]
    }
}

fn ik_jpsd(jacobian: Matrix, alpha: f64, initial: Matrix ,target: Matrix,  iteration: usize) -> Matrix {
    let res = Matrix::zero(jacobian.cols, 1);
    let mut jacobian_iter = jacobian;
    let mut initial_iter = initial;
    let mut angle_iter = Matrix {
        rows: 3,
        cols: 1,
        value: vec![
            vec![90.0*PI/180.0],
            vec![0.0],
            vec![0.0]
        ]
    };
    for i in 0..iteration {
            let mut parameter = Matrix::dot_mult(&jacobian_iter, &Matrix::mat_t(&jacobian_iter));
            parameter = Matrix::mat_inv(&mut parameter);
            Matrix::display(&parameter);
            /* 
            parameter = Matrix::dot_mult(&Matrix::mat_t(&mut jacobian_iter), &parameter);
            parameter = Matrix::dot_mult(&parameter, &error(&initial_iter, &target)); 
            parameter = Matrix::mat_scale(&parameter , alpha);
            println!("Iteration: {}", i+1);
            angle_iter = Matrix::mat_add(&parameter, &angle_iter);
            for i in 0..parameter.rows {
                if i == 0 {
                    print!("t{}: {}", i+1, rad_deg(angle_iter.value[i][0]));
                } else {
                    print!("   t{}: {}", i+1, rad_deg(angle_iter.value[i][0]));
                }
            }
            jacobian_iter = jacob(8.0, 5.0, 3.0, angle_iter.value[0][0], angle_iter.value[1][0], angle_iter.value[2][0]);
            initial_iter = xy(8.0, 3.0, 5.0, angle_iter.value[0][0], angle_iter.value[1][0], angle_iter.value[2][0]);
            let xy_current = initial_iter.clone(); println!();
            println!("Current Position: {},{}", xy_current.value[0][0], xy_current.value[1][0]); println!();
            
            */
        }
    return res
}

fn xy(l1: f64, l2: f64, l3: f64, t1: f64, t2: f64, t3: f64) -> Matrix {
    return Matrix {
        rows: 2,
        cols: 1,
        value: vec![
            vec![l1*cos(t1)+l2*cos(t1+t2)+l3*cos(t1+t2+t3)],
            vec![l1*sin(t1)+l2*sin(t1+t2)+l3*sin(t1+t2+t3)]
        ]
    }
}

fn cos(x: f64) -> f64 { return f64::cos(x); }
fn sin(x: f64) -> f64 { return f64::sin(x); }

fn error(initial: &Matrix, target: &Matrix) -> Matrix {
    return Matrix::mat_sub(target, initial);
}

fn rad_deg(x: f64) -> f64 {
    return x*180.0/PI
}

#[derive(Clone)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    value: Vec<Vec<f64>>
}

#[allow(unused)]
pub struct Rem {
    length: Vec<f64>,
    angle: Vec<f64>,
    target: Matrix
}

impl Matrix {
    pub fn zero(rows: usize, cols: usize) -> Self {
        return Matrix {
            rows: rows,
            cols: cols,
            value: vec![vec![0.0; cols]; rows]
        };
    }

    pub fn display(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{:?} ", self.value[i][j])
            }
            println!()
        }
    }

    pub fn mat_add(&self, others: &Matrix) -> Matrix {
        let mut res = Self::zero(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                res.value[i][j] = self.value[i][j] + others.value[i][j]
            }
        }
        return res
    }

    pub fn mat_sub(&self, others: &Matrix) -> Matrix {
        let mut res = Self::zero(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                res.value[i][j] = self.value[i][j] - others.value[i][j]
            }
        }
        return res
    }

    pub fn mat_scale(&self, scale: f64) -> Matrix {
        let mut res = Self::zero(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                res.value[i][j] = self.value[i][j] * scale
            }
        }
        return res
    }

    pub fn dot_mult(&self, others: &Matrix) -> Matrix {
        if self.cols != others.rows {
            panic!("Error of size Matrix")
        };
        let mut res = Self::zero(self.rows, others.cols);
        for i in 0..self.rows {
            for j in 0..others.cols {
                let mut sum: f64 = 0.0;
                for k in 0..self.cols {
                    sum += self.value[i][k] * others.value[k][j];
                }
                res.value[i][j] = sum;
            }
        }
        return res

    }

    pub fn mat_inv(&mut self) -> Matrix {
        let mut res = Self::zero(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                if self.value[i][j] < 0.000000001 {
                    self.value[i][j] = 0.0;
                } else {
                    res.value[i][j] = 1.0/self.value[i][j]
                }
            }
        }
        return res
    }

    pub fn mat_t(&self) -> Matrix {
        let mut res = Self::zero(self.cols, self.rows);
        for i in 0..res.rows {
            for j in 0..res.cols {
                res.value[i][j] = self.value[j][i]
            }
        }
        return res
    }

    pub fn mat_deg(&self) -> Matrix {
        let mut res = Self::zero(self.rows, self.cols);
        for i in 0..res.rows {
            for j in 0..res.cols {
            res.value[i][j] = self.value[i][j]*(180.0/PI)
            }
        }
        return res
    }

    pub fn mat_ide(size: usize) -> Matrix {
        let mut res = Self::zero(size,size);
        for i in 0..size {
            res.value[i][i] = 1.0;
        }
        return res
    }
}