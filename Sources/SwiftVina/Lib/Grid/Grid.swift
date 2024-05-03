//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

struct Grid {
    
    @usableFromInline
    var m_init: vec
    
    @usableFromInline
    var m_range: vec
    
    private let m_dim_fl_minus_1: vec
    
    var m_factor_inv: vec = vec(1,1,1)
    var m_data: Array3D<fl> = Array3D()
    
    private var m_factor: vec = vec(1,1,1)
    
    @usableFromInline
    var initialized: Bool {
        return m_data.dim0() > 0 && m_data.dim1() > 0 && m_data.dim2() > 0
    }

    init() {
        self.m_init = zero_vec
        self.m_range = vec(1,1,1)
        self.m_dim_fl_minus_1 = vec(-1,-1,-1)
    }
    
    init(gd: grid_dims) { 
        m_data.resize(gd[0].n_voxels + 1, gd[1].n_voxels + 1, gd[2].n_voxels + 1, fillValue: 0) // number of sample points == n_voxels + 1
        m_init  = vec(gd[0].begin,  gd[1].begin,  gd[2].begin)
        m_range = vec(gd[0].span(), gd[1].span(), gd[2].span())
        assert(m_range[0] > 0)
        assert(m_range[1] > 0)
        assert(m_range[2] > 0)
        m_dim_fl_minus_1 = vec(fl(m_data.dim0() - 1),
                               fl(m_data.dim1() - 1),
                               fl(m_data.dim2() - 1))
        
        // MARK: Test using SIMD
//        m_factor = m_dim_fl_minus_1 / m_range
//        m_factor_inv = recip(m_factor)
        for i in 0..<3 {
            m_factor[i] = m_dim_fl_minus_1[i] / m_range[i]
            m_factor_inv[i] = 1 / m_factor[i]
        }
        
    }

    @inlinable
    func index_to_argument(x: sz, y: sz, z: sz) -> vec {
        return vec(m_init[0] + m_factor_inv[0] * fl(x),
                   m_init[1] + m_factor_inv[1] * fl(y),
                   m_init[2] + m_factor_inv[2] * fl(z))
    }

    @inlinable
    func evaluate(location: vec, slope: fl, c: fl) -> fl {
        return evaluate_aux(location: location, slope: slope, v: c)
    }

    @inlinable
    func evaluate(location: vec, slope: fl, c: fl, deriv: inout vec) -> fl { // sets deriv
        return evaluate_aux_deriv(location: location, slope: slope, v: c, deriv: &deriv)
    }

    fileprivate func evaluate_aux(location: vec, slope: fl, v: fl) -> fl {
        var s: vec = (location - m_init) * m_factor
        var miss: vec = zero_vec
        var region: simd_long3 = simd_long3(0,0,0)
        var a: ContiguousArray<sz> = [0, 0, 0]
        
        for i in 0..<3 {
            let current_m_data_dim = m_data.dim(sz(i))
            let m_dim_fl_minus_1_i = m_dim_fl_minus_1[i]
            
            if s[i] < 0 {
                miss[i] = -s[i]
                region[i] = -1
                a[i] = 0
                s[i] = 0
            } else if s[i] >= m_dim_fl_minus_1_i {
                miss[i] = s[i] - m_dim_fl_minus_1_i
                region[i] = 1
                assert(current_m_data_dim >= 2)
                a[i] = current_m_data_dim - 2
                s[i] = 1
            } else {
                region[i] = 0 
                a[i] = sz(s[i])
                s[i] -= fl(a[i])
            }
            assert(s[i] >= 0)
            assert(s[i] <= 1)
            assert(a[i] >= 0)
            assert(a[i]+1 < current_m_data_dim)
        }
        
        let penalty: fl = slope * simd_dot(miss, m_factor_inv)
        assert(penalty > -epsilon_fl)

        let x0: sz = a[0]
        let y0: sz = a[1]
        let z0: sz = a[2]

        let x1: sz = x0+1
        let y1: sz = y0+1
        let z1: sz = z0+1

        let f000: fl = m_data[x0, y0, z0]
        let f100: fl = m_data[x1, y0, z0]
        let f010: fl = m_data[x0, y1, z0]
        let f110: fl = m_data[x1, y1, z0]
        let f001: fl = m_data[x0, y0, z1]
        let f101: fl = m_data[x1, y0, z1]
        let f011: fl = m_data[x0, y1, z1]
        let f111: fl = m_data[x1, y1, z1]

        let x: fl = s[0]
        let y: fl = s[1]
        let z: fl = s[2]

        let mx: fl = 1-x
        let my: fl = 1-y
        let mz: fl = 1-z

        var f = f000 * mx * my * mz
        f += f100 * x * my * mz
        f += f010 * mx * y * mz
        f += f110 * x * y * mz
        f += f001 * mx * my * z
        f += f101 * x * my * z
        f += f011 * mx * y * z
        f += f111 * x * y * z

        // Assuming SIMD3<Float> compatible types for vectors
//        let sVec = SIMD3<Double>(s[0], s[1], s[2])
//        let oneMinusS = SIMD3<Double>(1, 1, 1) - sVec
//
//        // Access the eight corner values of the cube in m_data
//        let corners = [
//            m_data[x0, y0, z0], m_data[x1, y0, z0], m_data[x0, y1, z0], m_data[x1, y1, z0],
//            m_data[x0, y0, z1], m_data[x1, y0, z1], m_data[x0, y1, z1], m_data[x1, y1, z1]
//        ]
//
//        // Create SIMD vectors from corner values
//        let f0 = SIMD4<Double>(corners[0], corners[1], corners[2], corners[3])
//        let f1 = SIMD4<Double>(corners[4], corners[5], corners[6], corners[7])
//
//        // Vectorized multiplication for interpolation
//        let mxmy_mz = oneMinusS.x * oneMinusS.y * oneMinusS.z
//        let mxmy_z  = oneMinusS.x * oneMinusS.y * sVec.z
//        let mx_y_mz = oneMinusS.x * sVec.y * oneMinusS.z
//        let mx_y_z  = oneMinusS.x * sVec.y * sVec.z
//        let x_my_mz = sVec.x * oneMinusS.y * oneMinusS.z
//        let x_my_z  = sVec.x * oneMinusS.y * sVec.z
//        let x_y_mz  = sVec.x * sVec.y * oneMinusS.z
//        let x_y_z   = sVec.x * sVec.y * sVec.z
//
//        // Calculate interpolated value
//        var interpolatedF = mxmy_mz * f0[0] + mxmy_z * f0[1] + mx_y_mz * f0[2] + mx_y_z * f0[3] +
//                            x_my_mz * f1[0] + x_my_z * f1[1] + x_y_mz * f1[2] + x_y_z * f1[3]
//
//
//        curl(e: &interpolatedF, v: v)
        curl(e: &f, v: v)
//        return interpolatedF + penalty
        return f + penalty
    }
    
    fileprivate func evaluate_aux_deriv(location: vec, slope: fl, v: fl, deriv: inout vec) -> fl {
        var s: vec = (location - m_init) * m_factor
        var miss: vec = zero_vec
        var region: simd_long3 = simd_long3(0,0,0)
        var a: ContiguousArray<sz> = [0, 0, 0]
        
        for i in 0..<3 {
            let current_m_data_dim = m_data.dim(sz(i))
            let m_dim_fl_minus_1_i = m_dim_fl_minus_1[i]
            if s[i] < 0 {
                miss[i] = -s[i]
                region[i] = -1
                a[i] = 0
                s[i] = 0
            } else if s[i] >= m_dim_fl_minus_1_i {
                miss[i] = s[i] - m_dim_fl_minus_1_i
                region[i] = 1
                assert(current_m_data_dim >= 2)
                a[i] = current_m_data_dim - 2
                s[i] = 1
            } else {
                region[i] = 0
                a[i] = sz(s[i])
                s[i] -= fl(a[i])
            }
            assert(s[i] >= 0)
            assert(s[i] <= 1)
            assert(a[i] >= 0)
            assert(a[i]+1 < current_m_data_dim)
        }
        
        let penalty: fl = slope * simd_dot(miss, m_factor_inv)
        assert(penalty > -epsilon_fl)

        let x0: sz = a[0]
        let y0: sz = a[1]
        let z0: sz = a[2]

        let x1: sz = x0+1
        let y1: sz = y0+1
        let z1: sz = z0+1

        let f000: fl = m_data[x0, y0, z0]
        let f100: fl = m_data[x1, y0, z0]
        let f010: fl = m_data[x0, y1, z0]
        let f110: fl = m_data[x1, y1, z0]
        let f001: fl = m_data[x0, y0, z1]
        let f101: fl = m_data[x1, y0, z1]
        let f011: fl = m_data[x0, y1, z1]
        let f111: fl = m_data[x1, y1, z1]

        let x: fl = s[0]
        let y: fl = s[1]
        let z: fl = s[2]

        let mx: fl = 1-x
        let my: fl = 1-y
        let mz: fl = 1-z

        var f  = f000 * mx * my * mz
            f += f100 * x * my * mz
            f += f010 * mx * y * mz
            f += f110 * x * y * mz
            f += f001 * mx * my * z
            f += f101 * x * my * z
            f += f011 * mx * y * z
            f += f111 * x * y * z

        let x_g = f000 * (-1)*my*mz + f100 * 1*my*mz + f010 * (-1)*y*mz + f110 * 1*y*mz + f001 * (-1)*my*z + f101 * 1*my*z + f011 * (-1)*y*z + f111 * 1*y*z
        let y_g = f000 * mx*(-1)*mz + f100 * x*(-1)*mz + f010 * mx*1*mz + f110 * x*1*mz + f001 * mx*(-1)*z + f101 * x*(-1)*z + f011 * mx*1*z + f111 * x*1*z
        let z_g = f000 * mx*my*(-1) + f100 * x*my*(-1) + f010 * mx*y*(-1) + f110 * x*y*(-1) + f001 * mx*my*1 + f101 * x*my*1 + f011 * mx*y*1 + f111 * x*y*1

        var gradient = vec(x_g, y_g, z_g)
        
        // Assuming SIMD3<Float> compatible types for vectors
//        let sVec = SIMD3<Double>(s[0], s[1], s[2])
//        let oneMinusS = SIMD3<Double>(1, 1, 1) - sVec
//
//        // Access the eight corner values of the cube in m_data
//        let corners = [
//            m_data[x0, y0, z0], m_data[x1, y0, z0], m_data[x0, y1, z0], m_data[x1, y1, z0],
//            m_data[x0, y0, z1], m_data[x1, y0, z1], m_data[x0, y1, z1], m_data[x1, y1, z1]
//        ]
//
//        // Create SIMD vectors from corner values
//        let f0 = SIMD4<Double>(corners[0], corners[1], corners[2], corners[3])
//        let f1 = SIMD4<Double>(corners[4], corners[5], corners[6], corners[7])
//
//        // Vectorized multiplication for interpolation
//        let mxmy_mz = oneMinusS.x * oneMinusS.y * oneMinusS.z
//        let mxmy_z  = oneMinusS.x * oneMinusS.y * sVec.z
//        let mx_y_mz = oneMinusS.x * sVec.y * oneMinusS.z
//        let mx_y_z  = oneMinusS.x * sVec.y * sVec.z
//        let x_my_mz = sVec.x * oneMinusS.y * oneMinusS.z
//        let x_my_z  = sVec.x * oneMinusS.y * sVec.z
//        let x_y_mz  = sVec.x * sVec.y * oneMinusS.z
//        let x_y_z   = sVec.x * sVec.y * sVec.z
//
//        // Calculate interpolated value
//        var interpolatedF = mxmy_mz * f0[0] + mxmy_z * f0[1] + mx_y_mz * f0[2] + mx_y_z * f0[3] +
//                            x_my_mz * f1[0] + x_my_z * f1[1] + x_y_mz * f1[2] + x_y_z * f1[3]
//        
//        // Calculate gradients for each dimension
//        let xGradient = SIMD4<Double>(-oneMinusS.y * oneMinusS.z, oneMinusS.y * oneMinusS.z, -sVec.y * oneMinusS.z, sVec.y * oneMinusS.z) * f0 +
//                        SIMD4<Double>(-oneMinusS.y * sVec.z, oneMinusS.y * sVec.z, -sVec.y * sVec.z, sVec.y * sVec.z) * f1
//
//        let yGradient = SIMD4<Double>(-oneMinusS.x * oneMinusS.z, -sVec.x * oneMinusS.z, oneMinusS.x * oneMinusS.z, sVec.x * oneMinusS.z) * f0 +
//                        SIMD4<Double>(-oneMinusS.x * sVec.z, -sVec.x * sVec.z, oneMinusS.x * sVec.z, sVec.x * sVec.z) * f1
//
//        let zGradient = SIMD4<Double>(-oneMinusS.x * oneMinusS.y, -sVec.x * oneMinusS.y, -oneMinusS.x * sVec.y, -sVec.x * sVec.y) * f0 +
//                        SIMD4<Double>(oneMinusS.x * oneMinusS.y, sVec.x * oneMinusS.y, oneMinusS.x * sVec.y, sVec.x * sVec.y) * f1
//
//        // Sum up the gradients to get the final gradient vector
//        
//        var gradient = SIMD3<Double>(reduce_add(xGradient), reduce_add(yGradient), reduce_add(zGradient))
//
//        curl(e: &interpolatedF, deriv: &gradient, v: v)
        
//        deriv = m_factor * (gradient * abs(region)) + slope * region
        
//        let gradient_everywhere: vec = vec(((region[0] == 0) ? gradient[0] : 0),
//                                           ((region[1] == 0) ? gradient[1] : 0),
//                                           ((region[2] == 0) ? gradient[2] : 0))
//
        curl(e: &f, deriv: &gradient, v: v)
        
        var gradient_everywhere: vec = zero_vec
        for i in 0..<3 {
            gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0)
            deriv[i] = m_factor[i] * gradient_everywhere[i] + slope * fl(region[i])
        }
        
//        deriv = m_factor * simd_select(gradient, zero_vec, region) + slope * vec(region)
        
//        return interpolatedF + penalty
        return f + penalty
    }
    
}
