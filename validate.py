'''
Validation and Benchmarking Script for Stats Functions Optimization

author: Qwen 3.5 397B (2026.11.15)
'''
import numpy as np
import time
from typing import Tuple, Dict, List

class CodeValidator:
    """Validates optimized code against original implementation."""
    
    def __init__(self, tolerance=1e-6):
        self.tolerance = tolerance
        self.results = {}
    
    def compare_arrays(self, name: str, old: np.ndarray, new: np.ndarray) -> bool:
        """Compare two arrays within tolerance."""
        if old.shape != new.shape:
            print(f"❌ {name}: Shape mismatch - Old: {old.shape}, New: {new.shape}")
            return False
        
        max_diff = np.max(np.abs(old - new))
        rel_diff = max_diff / (np.max(np.abs(old)) + 1e-10)
        
        if max_diff < self.tolerance or rel_diff < self.tolerance:
            print(f"✓ {name}: PASS (max_diff={max_diff:.2e})")
            return True
        else:
            print(f"❌ {name}: FAIL (max_diff={max_diff:.2e}, rel_diff={rel_diff:.2e})")
            return False
    
    def compare_scalars(self, name: str, old: float, new: float) -> bool:
        """Compare two scalar values within tolerance."""
        diff = abs(old - new)
        rel_diff = diff / (abs(old) + 1e-10)
        
        if diff < self.tolerance or rel_diff < self.tolerance:
            print(f"✓ {name}: PASS (diff={diff:.2e})")
            return True
        else:
            print(f"❌ {name}: FAIL (diff={diff:.2e})")
            return False


def run_validation_tests() -> Dict:
    """Run comprehensive validation tests."""
    validator = CodeValidator(tolerance=1e-4)  # Relaxed tolerance for FC method
    results = {'passed': 0, 'failed': 0, 'tests': []}
    
    print("=" * 70)
    print("VALIDATION TESTS: Original vs Optimized Code")
    print("=" * 70)
    
    # Test cases covering different regimes
    test_cases = [
        # (gross, bkg, p_conf, description)
        (5, 2, 0.68, "Low count detection"),
        (50, 30, 0.68, "Medium count detection"),
        (300, 200, 0.68, "High count regime"),
        (10, 15, 0.68, "Non-detection (bkg > gross)"),
        (100, 100, 0.95, "Medium count, 95% confidence"),
        (5, 2, 0.95, "Low count, 95% confidence"),
    ]
    
    for gross, bkg, p_conf, desc in test_cases:
        print(f"\nTest: {desc} (gross={gross}, bkg={bkg}, p_conf={p_conf})")
        print("-" * 70)
        
        # Test calc_source
        try:
            old_src = calc_source_old(gross, bkg)
            new_src = calc_source(gross, bkg)
            
            for i, name in enumerate(['src', 'src_lower', 'src_upper']):
                if validator.compare_scalars(f"calc_source.{name}", old_src[i], new_src[i]):
                    results['passed'] += 1
                else:
                    results['failed'] += 1
                results['tests'].append(f"calc_source.{name}_{desc}")
        except Exception as e:
            print(f"❌ calc_source failed: {e}")
            results['failed'] += 3
        
        # Test calc_conf
        try:
            old_conf = calc_conf_old(gross, bkg, p_conf)
            new_conf = calc_conf(gross, bkg, p_conf)
            
            for i, name in enumerate(['conf_lower', 'conf_upper']):
                if validator.compare_arrays(f"calc_conf.{name}", old_conf[i], new_conf[i]):
                    results['passed'] += 1
                else:
                    results['failed'] += 1
                results['tests'].append(f"calc_conf.{name}_{desc}")
        except Exception as e:
            print(f"❌ calc_conf failed: {e}")
            results['failed'] += 2
        
        # Test poiss_conf_lim
        try:
            old_lim = poiss_conf_lim_old(gross, np.array([0.1, 0.5, 0.9]))
            new_lim = poiss_conf_lim(gross, np.array([0.1, 0.5, 0.9]))
            
            if validator.compare_arrays(f"poiss_conf_lim", old_lim, new_lim):
                results['passed'] += 1
            else:
                results['failed'] += 1
            results['tests'].append(f"poiss_conf_lim_{desc}")
        except Exception as e:
            print(f"❌ poiss_conf_lim failed: {e}")
            results['failed'] += 1
    
    # Test feldman_cousins (most critical)
    print("\n" + "=" * 70)
    print("FELDMAN-COUSINS VALIDATION (Critical)")
    print("=" * 70)
    
    fc_cases = [
        (10, 5, 0.68),
        (20, 10, 0.68),
        (50, 30, 0.95),
    ]
    
    for gross, bkg, p_conf in fc_cases:
        try:
            old_fc = feldman_cousins_old(gross, bkg, p_conf)
            new_fc = feldman_cousins(gross, bkg, p_conf)
            
            for i, name in enumerate(['conf_lower', 'conf_upper']):
                if validator.compare_arrays(f"feldman_cousins.{name}", old_fc[i], new_fc[i]):
                    results['passed'] += 1
                else:
                    results['failed'] += 1
                results['tests'].append(f"feldman_cousins.{name}_g{gross}_b{bkg}")
        except Exception as e:
            print(f"❌ feldman_cousins failed: {e}")
            results['failed'] += 2
    
    print("\n" + "=" * 70)
    print(f"SUMMARY: {results['passed']} passed, {results['failed']} failed")
    print("=" * 70)
    
    return results


def run_benchmarks(n_iterations=100) -> Dict:
    """Benchmark performance improvements."""
    print("\n" + "=" * 70)
    print("PERFORMANCE BENCHMARKS")
    print("=" * 70)
    
    benchmarks = {}
    
    # Benchmark calc_source
    test_params = [(10, 5), (50, 30), (300, 200)]
    
    for gross, bkg in test_params:
        # Old version
        start = time.perf_counter()
        for _ in range(n_iterations):
            calc_source_old(gross, bkg)
        old_time = (time.perf_counter() - start) / n_iterations
        
        # New version
        start = time.perf_counter()
        for _ in range(n_iterations):
            calc_source(gross, bkg)
        new_time = (time.perf_counter() - start) / n_iterations
        
        speedup = old_time / new_time if new_time > 0 else float('inf')
        
        print(f"\ncalc_source (gross={gross}, bkg={bkg}):")
        print(f"  Old: {old_time*1000:.3f} ms")
        print(f"  New: {new_time*1000:.3f} ms")
        print(f"  Speedup: {speedup:.2f}x")
        
        benchmarks[f"calc_source_g{gross}_b{bkg}"] = {
            'old_ms': old_time * 1000,
            'new_ms': new_time * 1000,
            'speedup': speedup
        }
    
    # Benchmark feldman_cousins (most significant improvement)
    print("\n" + "-" * 70)
    print("feldman_cousins (Most Critical Optimization):")
    print("-" * 70)
    
    for gross, bkg in [(20, 10), (50, 30)]:
        # Old version
        start = time.perf_counter()
        for _ in range(10):  # Fewer iterations for FC
            feldman_cousins_old(gross, bkg, 0.68)
        old_time = (time.perf_counter() - start) / 10
        
        # New version
        start = time.perf_counter()
        for _ in range(10):
            feldman_cousins(gross, bkg, 0.68)
        new_time = (time.perf_counter() - start) / 10
        
        speedup = old_time / new_time if new_time > 0 else float('inf')
        
        print(f"\nfeldman_cousins (gross={gross}, bkg={bkg}):")
        print(f"  Old: {old_time*1000:.3f} ms")
        print(f"  New: {new_time*1000:.3f} ms")
        print(f"  Speedup: {speedup:.2f}x")
        
        benchmarks[f"feldman_cousins_g{gross}_b{bkg}"] = {
            'old_ms': old_time * 1000,
            'new_ms': new_time * 1000,
            'speedup': speedup
        }
    
    return benchmarks


if __name__ == "__main__":
    # Run validation
    validation_results = run_validation_tests()
    
    # Run benchmarks
    benchmark_results = run_benchmarks()
    
    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Validation: {validation_results['passed']}/{validation_results['passed'] + validation_results['failed']} tests passed")
    print(f"Average speedup: {np.mean([b['speedup'] for b in benchmark_results.values()]):.2f}x")
    print("=" * 70)
