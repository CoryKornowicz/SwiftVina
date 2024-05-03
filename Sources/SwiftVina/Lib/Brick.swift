//
//  Brick.swift
//
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

// inline fl closest_between(fl begin, fl end, fl x) {
// 	assert(begin <= end);
// 	if(x <= begin) return begin;
// 	else if(x >= end) return end;
// 	return x;
// }

// inline vec brick_closest(const vec& begin, const vec& end, const vec& v) {
// 	vec tmp;
// 	VINA_FOR_IN(i, tmp)
// 		tmp[i] = closest_between(begin[i], end[i], v[i]);
// 	return tmp;
// }

// inline fl brick_distance_sqr(const vec& begin, const vec& end, const vec& v) {
// 	vec closest; closest = brick_closest(begin, end, v);
// 	return vec_distance_sqr(closest, v);
// }

@inlinable
func closest_between(_ begin: fl, _ end: fl, _ x: fl) -> fl {
    assert(begin <= end)
    if x <= begin {
        return begin
    }
    else if x >= end {
        return end
    }
    return x
}

@inlinable
func brick_closest(_ begin: vec, _ end: vec, _ v: vec) -> vec {
    var tmp: vec = zero_vec
    for i in 0..<3 {
        tmp[i] = closest_between(begin[i], end[i], v[i])
    }
    return tmp
}

//@inlinable
//func brick_closest(_ begin: vec, _ end: vec, _ v: vec) -> vec {
//    let minVec = min(begin, end)
//    let maxVec = max(begin, end)
//    return max(minVec, min(maxVec, v))
//}

@inlinable
func brick_distance_sqr(begin: vec, end: vec, v: vec) -> fl {
    let closest = brick_closest(begin, end, v)
    return vec_distance_sqr(closest, v)
}
